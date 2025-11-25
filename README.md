# Mixed-Strategy Primary and Backup Controller Optimization

### Instructions

Ensure that you have both [AMPL](https://ampl.com/) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) installed on your system.

Create a virtual environment

```bash
uv venv --python 3.13 $(pwd)
source .venv/bin/activate
uv sync
```

Create a `.env` file and add the environment variable for your AMPL license UUID.

```
AMPL_LICENSE_UUID=xxxxxxxxxxx
```

Enjoy

### Creating a new model

To create a new model, we use the AMPL language.

The class must be:

```python
class NewNetworkModel(MathematicalModel):
    def __init__(self, network, **params):
        self._model_file = 'models/my_new_network_model.mod'
        super(NewNetworkmodel, self).__init__(self._model_file)
        self._name = 'New Network Model'
        self.network = network
    def report(self):
        pass
        # Should return the solved answers
    def load_data(self, **params):
        pass
        # Load data into ampl.
```
