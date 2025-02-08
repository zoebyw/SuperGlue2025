from flask import Flask, render_template, jsonify
from flask_cors import CORS
from config import config
from load_data import DataLoader
app = Flask(__name__, template_folder='templates')
app.config.from_object(config['default'])
CORS(app)
data_loader = DataLoader()
@app.route('/')
# def index():
#     return render_template('index.html')
@app.route('/compounds', methods=['GET'])
def get_compounds():
    # try:
    #     compounds = data_loader.read_compounds()
    #     return jsonify({
    #         'status': 'success',
    #         'compounds': compounds
    #     })
    # except Exception as e:
    #     return jsonify({
    #         'status': 'error',
    #         'message': str(e)
    #     }), 400
    try:
        compounds = data_loader.read_compounds()
        # print(compounds)
        return render_template('index.html', compounds=compounds)
    except Exception as e:
        return f"Error: {str(e)}", 500
if __name__ == '__main__':
    app.run(debug=True)