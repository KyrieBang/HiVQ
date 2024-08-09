from flask import Flask,send_file,render_template


app = Flask(__name__)

@app.route("/HiVQ")
def returnHtml():
    return render_template("HiVQ.html")


if __name__ == '__main__':
    app.run(host='0.0.0.0', port='8888', debug=True)