

import mysql.connector
from hevelius import config

def connect():
    try:
        cnx = mysql.connector.connect(user=config.USER, password=config.PASSWORD, database=config.DBNAME, host=config.HOST, port=config.PORT)
    except BaseException as e:
        print("ERROR: Failed to connect to DB: user=%s, database=%s, host=%s, port=%d: %s" % (config.USER, config.DBNAME, config.HOST, config.PORT, e))
        raise
    return cnx

def run_query(cnx, query):
    cursor = cnx.cursor() # cursor(dictionary=True) or cursor(named_tuple=True)
    cursor.execute(query)
    try:
        result = cursor.fetchall()
    except mysql.connector.Error as err:
        print("#### Something went wrong: {}".format(err))
#    except:
#        #result = None # If this is an update or delete query.
#        cnx.commit()
    cnx.commit()
    cursor.close()
    return result

def version_get(cnx):
    query = 'SELECT * from schema_version'
    cursor = cnx.cursor()

    try:
        cursor.execute(query)
    except mysql.connector.errors.ProgrammingError:
        # Table doesn't exist, return 0
        return 0

    for i in cursor:
        ver = i[0]

    cursor.close()
    return ver

def stats_print(cnx):

    stats = [
        ("true", "all tasks"),
        ("imagename is not null", "with images"),
        ("he_solved = 1",            "solved"),
        ("he_solved is null",        "not solved (not attempted)"),
        ("he_solved = 0",            "not solved (attempted, but failed)"),
        ("he_fwhm is not null", "with quality measured (FWHM != null)"),
        ("eccentricity is not null", "with quality measured (eccen != null)"),

    ]


    for cond, descr in stats:
        q = "SELECT count(*) FROM tasks WHERE %s" % cond

        result = run_query(cnx, q)[0][0]

        print("Tasks %40s: %d" % (descr, result))

    # Get overall tasks counters
    # tasks_cnt = run_query(cnx, 'SELECT count(*) from tasks')[0][0]
    # files_cnt = run_query(cnx, 'SELECT count(*) from tasks where imagename is not null')[0][0]

    # fwhm_cnt = run_query(cnx, 'select count(*) from tasks WHERE fwhm is not null')
    # eccentricity_cnt = run_query(cnx, 'select count(*) from tasks WHERE eccentricity is not null')

    # print("There are %d tasks, %d files, %d have FWHM, %d have eccentricity." % stats)
    # print("Missing: %d files miss FWHM, %d files miss eccentricity." % (stats[1] - stats[2], stats[1] - stats[3]))


    # return tasks_cnt[0][0], files_cnt[0][0], fwhm_cnt[0][0], eccentricity_cnt[0][0]

def stats_by_state(cnx):
    # Get tasks list by status
    hist = run_query(cnx, 'SELECT state, name, count(*) from tasks, states where tasks.state = states.id group by state')
    res = []

    for row in hist:
        res.append( (row[0], row[1], row[2]))

    return res

def stats_by_user(cnx, state = 6):
    if state is None:
        cond = ""
    else:
        cond = "AND state = %d" % state

    q = "SELECT login, tasks.user_id, count(*) from tasks, users where tasks.user_id = users.user_id %s group by tasks.user_id order by login" % cond

    tasks_per_user = run_query(cnx, q)

    res = []
    for row in tasks_per_user:
        res.append( (row[0], row[1], row[2]) )
    return res

def task_get(cnx, id):
    q = "SELECT task_id, state, user_id, imagename, object, descr, comment, ra, decl, exposure, filter, binning, guiding, fwhm, eccentricity FROM tasks WHERE task_id = %d" % id

    t = run_query(cnx, q)[0]

    x = {}
    x["id"] = t[0]
    x["state"] = t[1]
    x["user_id"] = t[2]
    x["file"] = t[3]
    x["object"] = t[4]
    x["descr"] = t[5]
    x["comment"] = t[6]
    x["ra"] = t[7]
    x["decl"] = t[8]
    x["exposure"] = t[9]
    x["filter"] = t[10]
    x["binning"] = t[11]
    x["guiding"] = t[12]
    x["fwhm"] = t[13]
    x["eccentricity"] = t[14]

    return x

def tasks_get_filter(cnx, criteria):
    q = f"SELECT state,task_id, imagename, object, he_solved_ra, he_solved_dec, exposure, filter, binning, fwhm, eccentricity FROM tasks WHERE {criteria}"

    tasks = run_query(cnx, q)

    print(f"Selected {len(tasks)} task(s)")

    return tasks

def task_update(cnx, id, fwhm = None, eccentricity = None):
    upd = ""
    if fwhm is not None:
        upd = "fwhm = %f" % fwhm
    if eccentricity is not None:
        if len(upd):
            upd += ", "
        upd += "eccentricity = %f" % eccentricity

    if not len(upd):
        print("Nothing to update in task %d, aborting" % id)

    q = "UPDATE tasks SET %s WHERE task_id=%d" % (upd, id)

    print("Updating task %d: query=[%s]" % (id, q))

    run_query(cnx, q)