"""
Flask application that provides a REST API to the Hevelius backend.
"""

from flask import Flask, render_template, request
from flask_cors import CORS, cross_origin
import pandas as pd
import json
import plotly
import plotly.express as px

from hevelius import cmd_stats, db

# By default, Flask searches for templates in the templates/ dir.
# Other params: debug=True, port=8080
app = Flask(__name__)
CORS(app)


@app.route('/')
def root():
    """Just a stub API homepage."""
    return "Home🏠 EE"


@app.route('/histo')
def histogram():
    """Generates 2D diagram of observation density. Returns a HTML page with
    embedded plotly image."""

    # example data input
    # df = pd.DataFrame({
    #     'Fruit': ['Apples', 'Oranges', 'Bananas', 'Apples', 'Oranges',
    #               'Bananas'],
    #     'Amount': [4, 1, 2, 2, 4, 5],
    #     'City': ['SF', 'SF', 'SF', 'Montreal', 'Montreal', 'Montreal']
    # })
    # fig = px.bar(df, x='Fruit', y='Amount', color='City', barmode='group')

    histo = cmd_stats.histogram({})

    x_labels = [cmd_stats.deg2rah(a) for a in range(0, 360, 1)]
    y_labels = [str(y) for y in range(90, -90, -1)]

    labels = dict(x="Right Ascension (h:m/deg)",
                  y="Declination (deg)", color="# of frames")

    pandas = pd.DataFrame(histo)
    fig = px.imshow(pandas, labels=labels, y=y_labels, x=x_labels)

    fig['layout']['yaxis']['autorange'] = "reversed"

    graph_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return render_template('histogram.html', graphJSON=graph_json)


@app.route('/api/tasks', methods=['POST', 'GET'])
def tasks():
    """
    Flask method used when a list of tasks is required.
    """

    user_id = get_param(request, 'user_id')
    limit = get_param(request, 'limit')

    q = "SELECT task_id, tasks.user_id, aavso_id, object, ra, decl, " \
        "exposure, descr, filter, binning, guiding, dither, " \
        "defocus, calibrate, solve, vphot, other_cmd, " \
        "min_alt, moon_distance, skip_before, skip_after, " \
        "min_interval, comment, state, imagename, " \
        "created, activated, performed, max_moon_phase, " \
        "max_sun_alt, auto_center, calibrated, solved, " \
        "sent FROM tasks, users WHERE tasks.user_id = users.user_id"
    if user_id is not None:
        q = q + f" AND tasks.user_id={user_id}"

    q = q + " ORDER by task_id DESC"

    if limit is not None:
        q = q + f" LIMIT {limit}"

    cnx = db.connect()
    tasks_list = db.run_query(cnx, q)
    cnx.close()

    return tasks_list


def sanitize(x: str) -> str:
    """
    Sanitizes x input (removes backslashes)
    """
    x = str(x).replace('\'', '')  # apostrophes are bad
    x = x.replace(';', '')  # commas also
    x = x.replace('\\', '')  # and backslashes
    return x


def get_param(req, field) -> str:
    json_html_request = req.get_json()
    x = json_html_request.get(field)
    if x:
        return sanitize(x)
    return x


@app.route('/api/login', methods=['POST'])
@cross_origin()
def login():
    """

    Example request:    {"username":"tomek","password":"1264a4e31234abcdef"}

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
    """

    if not request.is_json:
        # The input is totally messed up.
        return {'status': False, 'msg': 'Need JSON input'}

    user = get_param(request, 'username')
    md5pass = get_param(request, 'password')

    if user is None:
        return {'status': False, 'msg': 'Username not provided'}
    if md5pass is None:
        return {'status': False, 'msg': 'Password not provided'}

    q = f"""SELECT user_id, pass_d, login, firstname, lastname, share, phone, email, permissions,
            aavso_id, ftp_login, ftp_pass FROM users WHERE login='{user}'"""

    cnx = db.connect()
    db_resp = db.run_query(cnx, q)
    cnx.close()

    if db_resp is None or not len(db_resp):
        print(f"Login: No such username ({user}")
        return {'status': False, 'msg': 'Invalid credentials'}

    user_id, pass_db, _, firstname, lastname, share, phone, email, permissions, aavso_id, \
        ftp_login, ftp_pass = db_resp[0]

    if md5pass.lower() != pass_db.lower():
        print(f"Login ({user} exists, invalid pwd: expected {pass_db.lower()}, got {md5pass.lower()}")
        # Password's MD5 did not match
        return {'status': False, 'msg': 'Invalid credentials'}

    print(f"User {user} provided valid password, logged in.")
    return {'status': True,
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
            'msg': 'Welcome'}
