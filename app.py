from flask import Flask, render_template, jsonify

app = Flask(__name__, static_folder='static', template_folder='templates')

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')
if __name__ == '__main__':
    app.run(debug=True)
