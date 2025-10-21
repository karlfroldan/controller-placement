module SNDlibParser

export SndLibNode, SndLibEdge, SndLibGraph, parse_sndlib_graph

struct SndLibNode
    id  :: String
    lon :: Float64
    lat :: Float64
end

struct SndLibEdge
    id     :: String
    source :: String
    target :: String
end

struct SndLibGraph
    nodes :: Vector{SndLibNode}
    edges :: Vector{SndLibEdge}
end

"""
    parse_sndlib_graph(text::AbstractString) -> SndLibGraph

No validation. Parses only NODES and LINKS sections.
"""
function parse_sndlib_graph(text::AbstractString)
    nodes_block = _extract_block(text, "NODES")
    links_block = _extract_block(text, "LINKS")

    nodes = _parse_nodes(nodes_block)
    edges = _parse_links(links_block)

    SndLibGraph(nodes, edges)
end

# --- internals ---

# Robust balanced-paren extractor:
# finds the first '(' after the header token and returns its balanced contents.
function _extract_block(text::AbstractString, header::AbstractString)
    hstart = findfirst(Regex("\\b" * header * "\\b"), text)
    hstart === nothing && return ""
    open_idx = findnext('(', text, first(hstart))
    open_idx === nothing && return ""

    depth = 0
    content_start = 0
    i = open_idx
    while i <= lastindex(text)
        c = text[i]
        if c == '('
            depth += 1
            if depth == 1
                content_start = nextind(text, i)  # first char after '('
            end
        elseif c == ')'
            depth -= 1
            if depth == 0
                return text[content_start:prevind(text, i)]
            end
        end
        i = nextind(text, i)
    end
    return ""  # unbalanced -> empty
end

function _parse_nodes(block::AbstractString)
    nodes = SndLibNode[]
    for line in eachline(IOBuffer(block))
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        # Example: "Amsterdam ( 4.90 52.35 )"
        m = match(r"^([A-Za-z]+)\s*\(\s*([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s*\)$", s)
        m === nothing && continue
        city = m.captures[1]
        lon  = parse(Float64, m.captures[2])
        lat  = parse(Float64, m.captures[3])
        push!(nodes, SndLibNode(city, lon, lat))
    end
    nodes
end

function _parse_links(block::AbstractString)
    edges = SndLibEdge[]
    for line in eachline(IOBuffer(block))
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        
        m = match(r"^([A-Za-z0-9]+)\s*\(\s*([A-Za-z]+)\s+([A-Za-z]+)\s*\)", s)
        m === nothing && continue
        link_id = m.captures[1]
        src     = m.captures[2]
        dst     = m.captures[3]
        push!(edges, SndLibEdge(link_id, src, dst))
    end
    edges
end

end
