from flask import Flask, render_template
from flask_cors import CORS
from config import config
app = Flask(__name__, template_folder='templates')
app.config.from_object(config['default'])
CORS(app)

@app.route('/')
def index():
    return render_template('index.html')
if __name__ == '__main__':
    app.run(debug=True)