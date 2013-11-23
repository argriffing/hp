"""
Write some json output given an input file with a bunch of search results.
"""

import json

def gen_paragraphs(fin):
    paragraph = []
    for line in fin:
        if not line.strip():
            if paragraph:
                yield paragraph
                paragraph = []
        else:
            paragraph.append(line)
    if paragraph:
        yield paragraph

def write_json(fin, fout):
    """
    The json output will be over-processed.

    This over-processing is because I am more fluent in python
    than in javascript, so I will do things like follow references
    in python so that I can avoid having to follow references in javascript.

    """
    hp_to_value = {'h' : 1, 'p' : 0}
    direction_to_delta = {
            'l': (-1, 0),
            'r': (1, 0),
            'u': (0, -1),
            'd': (0, 1),
            }
    cells = []
    for paragraph in gen_paragraphs(fin):
        
        # Read the cell info from the paragraph.
        d = dict(line.split() for line in paragraph)
        n = int(d['n'])
        k = int(d['k'])
        score = int(d['score'])
        conformation = d['conformation']
        hp_assignment = d['hp']

        # Initialize the cell dictionary.
        json_nodes = []
        json_primary_edges = []
        json_spatial_edges = []
        cell = {
                'n': n,
                'k': k,
                'score': score,
                'nodes': json_nodes,
                'primaryEdges': json_primary_edges,
                'spatialEdges': json_spatial_edges,
                }
        cells.append(cell)

        # Follow the steps to construct the conformation.
        origin = (0, 0)
        points = [origin]
        spatial_edges = []
        primary_edges = []
        occupants = {origin : 0}
        for step_index, direction in enumerate(conformation):
            dx, dy = direction_to_delta[direction]
            pprev = points[-1]
            xprev, yprev = pprev
            xcurr, ycurr = xprev + dx, yprev + dy
            pcurr = (xcurr, ycurr)
            points.append(pcurr)
            primary_edges.append((pprev, pcurr))
            occupants[pcurr] = step_index + 1
            for dx, dy in direction_to_delta.values():
                nx = xcurr + dx
                ny = ycurr + dy
                pneighbor = (nx, ny)
                neighbor = occupants.get(pneighbor, None)
                if neighbor is not None and neighbor < step_index:
                    if hp_assignment[neighbor] == 'h':
                        if hp_assignment[step_index + 1] == 'h':
                            spatial_edges.append((pcurr, pneighbor))

        # Fill the json arrays.
        for (x, y), hp in zip(points, hp_assignment):
            json_nodes.append({
                'x': x,
                'y': y,
                'value': hp_to_value[hp],
                })
        for ((x1, y1), (x2, y2)) in primary_edges:
            json_primary_edges.append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
        for ((x1, y1), (x2, y2)) in spatial_edges:
            json_spatial_edges.append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
    
    # Write the cell array in json format.
    fout.write(json.dumps(cells, indent=2))


def main():
    json_filename = 'sol.json'
    solutions_filename = 'nksolutions.out'
    with open(solutions_filename) as fin:
        with open(json_filename, 'w') as fout:
            write_json(fin, fout)

if __name__ == '__main__':
    main()

