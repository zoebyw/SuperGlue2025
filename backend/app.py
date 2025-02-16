from flask import Flask
from flask_cors import CORS
from molecule_annotate import get_compounds, get_compound
from file_upload import upload_file
from molecule_convert import convert_molecule
app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

@app.route('/api/upload', methods=['POST'])
def handle_upload():
    return upload_file()
@app.route('/api/convert_molecule', methods=['POST'])
def handle_convert_molecule():
    return convert_molecule()
@app.route('/api/compounds', methods=['GET'])
def handle_get_compounds():
    return get_compounds()

# @app.route('/api/compound/<cmpd_id>', methods=['GET'])
# def handle_get_compound(cmpd_id):
#     return get_compound(cmpd_id)

if __name__ == '__main__':
    app.run(debug=True, port=5001)