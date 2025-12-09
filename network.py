from itertools import product
import numpy as np
import networkx as nx
from networkx import Graph

import matplotlib.pyplot as plt

network_files = {
    'dognet': 'dognet.dat',
    'cost266': 'net-cost266-l.dat',
    'coronet-conus': 'coronet_conus-l.dat',
}

class Network:
    def __init__(self, network = None):
        self.g = Graph()

        if network is not None and isinstance(network, str):
            filename = network_files[network]
            self.load_from_file(f'networks/{filename}')

        
        if network:
            n = self.order()
            self.delays = np.zeros((n, n))
            self._calculate_delays()

    def __getattr__(self, name):
        return getattr(self.g, name)

    def __str__(self):
        return str(self.g)

    def attack(self, attack_ids):
        """Attack the network given the attack ids. Returns a list of the connected components"""

        g_nodes = set(self.g)

        if isinstance(attack_ids, int):
            attack_ids = [attack_ids]

        induce = lambda nodes: nx.induced_subgraph(self.g, nodes)

        remaining_nodes = g_nodes.difference(set(attack_ids))

        H = induce(remaining_nodes)
        # This recreates the remaining component subgraphs
        connected_components = map(lambda s: induce(sorted(list(s))), nx.connected_components(H))
        return list(connected_components)

    def surviving_nodes(self, controller_ids, attack_ids, backup_controllers=None):
        connected_components = self.attack(attack_ids)
        survives = []

        # A component survives iff it has a controller
        controller_ids = set(controller_ids)
        if backup_controllers:
            backup_controllers = set(backup_controllers)
            controller_ids = controller_ids.union(backup_controllers)
            
        # backup_controllers = set(backup_controllers)
        for c in connected_components:
            c_nodes = set(c)
            if c_nodes.intersection(controller_ids):
                survives.extend(c_nodes)
        return survives

    def node_location(self, node_id):
        """Return the location of node `node_id` as a tuple (longitude, latitude)"""
        if self.has_node(node_id):
            return (self.g.nodes[node_id]['loc_x'], self.g.nodes[node_id]['loc_y'])
        raise Exception(f'node {node_id} does not exist in the graph')

    def _calculate_delays(self):
        """Calculate the network delay using dijkstra"""
        for v in self.g:
            # The first index of the tuple are the distances to go from one node to another.
            path_lengths = nx.single_source_dijkstra(self.g, v, weight='weight')[0]
            for w, dist in path_lengths.items():
                self.delays[v, w] = dist

    def draw(self, node_id='name'):
        pos = {}
        labels = {}
        for ix, n in enumerate(self.g.nodes):
            p = (self.g.nodes[n]['loc_y'], self.g.nodes[n]['loc_x'])
            pos[n] = p
            if node_id == 'name':
                labels[n] = self.g.nodes[n]['label']
            else:
                labels[n] = ix
        nx.draw(self.g, pos=pos, with_labels=True, labels=labels)
        plt.show()

    def load_from_file(self, filename):
        parsed_data = {
            'nodes': [],
            'edges': [],
        }

        on_node = False
        on_edge = False

        with open(filename, 'r') as f:
            raw_data = f.read()

        for line in raw_data.splitlines():
            if line.startswith('param: sa_Nodes:'):
                on_node = True
                continue

            if line.startswith('param: sa_Links:'):
                on_edge = True
                continue

            if line.startswith(';'):
                if on_node:
                    on_node = False
                else:
                    on_edge = False

            if on_node:
                parts = line.split()

                if int(parts[0]) <= 0:
                    raise Exception('Node ID should be 1-indexed')

                parsed_data['nodes'].append({
                    'id': int(parts[0]),
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'label': str(parts[3]),
                })

                continue

            if on_edge:
                parts = line.split()

                parsed_data['edges'].append({
                    'id': int(parts[0]),
                    'src': int(parts[1]),
                    'dst': int(parts[2]),
                    'w': float(parts[3]),
                })

                continue
        if on_edge and on_node:
            raise Exception('Malformed Data File')

        for n in parsed_data['nodes']:
            self.g.add_node(n['id'] - 1, loc_x=n['x'], loc_y=n['y'], label=n['label'])
        
        for e in parsed_data['edges']:
            self.g.add_edge(e['src'] - 1, e['dst'] - 1, weight=e['w'])
        
if __name__ == '__main__':
    n = Network('networks/dognet.dat')
    ax = n.surviving_nodes([1, 5], [5])
    print(ax)
