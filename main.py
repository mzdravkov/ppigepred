from flask import Flask, request
import json

import ppigepred.ui


if __name__ == '__main__':
    node_data, edges = ppigepred.ui.UI.run()

    app = Flask(__name__)

    @app.route('/')
    def visualize():
        f = open('ppigepred/visualize.html')
        return f.read()

    @app.route('/graphjs')
    def graphjs():
        f = open('ppigepred/graph.js')
        return f.read()

    @app.route('/subgraph')
    def serve_graph():
        data = {
            'nodes': node_data,
            'edges': edges,
        }
        return json.dumps(data)

    print('Visualization available at http://localhost:5000/')
    app.run(debug=True, use_reloader=False)
