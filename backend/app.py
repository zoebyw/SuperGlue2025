from flask import Flask,send_from_directory,jsonify,request
from flask_cors import CORS
from molecule_annotate import get_compounds, get_compound
from file_upload import upload_file
from molecule_convert import convert_molecule
from molecule_visualize import MoleculeVisualizer
import os
app = Flask(__name__)
CORS(app,resources={
    r"/api/*":{"origins":"*"},
    r"/data/*": {"origins": "*"},
    r"/get_molecule_image/*": {"origins": "*"}
    })
  

visualizer = MoleculeVisualizer(data_dir='data')
UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'data')
@app.route('/api/upload', methods=['POST'])
def handle_upload():
    return upload_file()

@app.route('/data/<filename>')
def serve_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)
@app.route('/api/convert_molecule', methods=['POST'])
def handle_convert_molecule():
    return convert_molecule()
@app.route('/api/compounds', methods=['GET'])
def handle_get_compounds():
    return get_compounds()

@app.route('/get_molecule_image/<cmpd_id>', methods=['GET'])
def handle_visualize(cmpd_id):
    try:
        filename = request.args.get('filename')
        if not filename:
            return jsonify({
                'success': False,
                'error': 'Missing filename'
            }), 400
        result = visualizer.process_request(cmpd_id, filename)
        return jsonify(result)
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

# @app.route('/api/compound/<cmpd_id>', methods=['GET'])
# def handle_get_compound(cmpd_id):
#     return get_compound(cmpd_id)

if __name__ == '__main__':
    app.run(debug=True, port=5001)