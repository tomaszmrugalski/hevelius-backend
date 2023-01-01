from flask import Flask, render_template, request
from flask_cors import CORS, cross_origin
import pandas as pd
import json
import plotly
import plotly.express as px

from hevelius import db_mysql as db

# By default, Flask searches for templates in the templates/ dir.
# Other params: debug=True, port=8080
app = Flask(__name__)
CORS(app)

@app.route('/')
def root():
    return "Home🏠 EE"

@app.route('/histo')
def histogram():
   df = pd.DataFrame({
      'Fruit': ['Apples', 'Oranges', 'Bananas', 'Apples', 'Oranges',
      'Bananas'],
      'Amount': [4, 1, 2, 2, 4, 5],
      'City': ['SF', 'SF', 'SF', 'Montreal', 'Montreal', 'Montreal']
   })

   fig = px.bar(df, x='Fruit', y='Amount', color='City', barmode='group')
   graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
   return render_template('histogram.html', graphJSON=graphJSON)

@app.route('/api/tasks')
def tasks():

    cnx = db.connect()
    tasks = db.tasks_get_filter(cnx, "imagename is not null AND he_solved_ra is not null AND state = 6 LIMIT 10")
    cnx.close()

    t = [ {
        "task_id": 123,
        "ra": 12.34,
        "decl": 45.67,
        "descr": "some object"
    }]

    return tasks


def sanitize(x: str) -> str:
    x = x.replace('\'','') # apostrophes are bad
    x = x.replace(';','') # commas also
    x = x.replace('\\','') # and backslashes
    return x

def get_param(request, field) -> str:
    json = request.get_json()
    x = json.get(field)
    if x:
        return sanitize(x)
    return x

@app.route('/api/login', methods = ['POST'])
@cross_origin()
def login():

    """_summary_

    Example request:    {"username":"tomek","password":"64a4e327e97c1e7926f9240edb123456"}

    Example response:
    {
        "aavso_id":"MTOA",
        "email":"some.mail@example.org",
        "firstname":"Tomek",
        "ftp_login":"MTOA",
        "ftp_pass":"xxxxx",
        "lastname":"M.",
        "msg":"Welcome",
        "permissions":1,
        "phone":"",
        "share":0.0,
        "status":true,
        "user_id":3
    }

    :return: _description_
    :rtype: _type_
    """

    if not request.is_json:
        return { 'status': False, 'msg': 'Need JSON input' }

    user = get_param(request, 'username')
    md5pass = get_param(request, 'password')

    if user is None:
        return { 'status': False, 'msg': 'Username not provided'}
    if md5pass is None:
        return { 'status': False, 'msg': 'Password not provided'}

    q = f"SELECT user_id, pass_d, login, firstname, lastname, share, phone, email, permissions, aavso_id, ftp_login, ftp_pass FROM users WHERE login='{user}'"

    cnx = db.connect()
    db_resp = db.run_query(cnx, q)
    cnx.close()

    if db_resp is None or not len(db_resp):
        print("#### No username")
        return { 'status': False, 'msg': 'Invalid credentials'} # No such username

    user_id, pass_db, login, firstname, lastname, share, phone, email, permissions, aavso_id, ftp_login, ftp_pass = db_resp[0]
    if (md5pass.lower() != pass_db.lower()):
        print(f"#### Login exists, but invalid password: expected {pass_db.lower()}, got {md5pass.lower()}")
        return { 'status': False, 'msg': 'Invalid credentials'} # Password's MD5 did not match

    return { 'status': True,
             'user_id': user_id,
             'firstname': firstname,
             'lastname': lastname,
             'share': share,
             'phone': phone,
             'email': email,
             'permissions': permissions,
             'aavso_id': aavso_id,
             'ftp_login': ftp_login,
             'ftp_pass': ftp_pass,
             'msg': 'Welcome' }