"""

This module is responsible for handling files from the iteleskop.org project.

"""

import re
import os


def parse_iteleskop_filename(fname: str) -> dict:
    """
    Parses name of the file that uses iteleskop naming convention:

    SFDB_(date)_J(task_id)_AVSO_FILTER_BIN_EXP_OBJECT.fit

    :param fname: filename of the file
    :return: dictionary with parsed parameters (or empty dict)
    """

    # Step 1: get rid of the possible leading dir name.
    dir, filename = os.path.split(fname)  # get rid of the directory
    filename, _ = os.path.splitext(filename)  # and the extension

    # pylint: disable=line-too-long
    #        (flags           ) (date                                                  )  (job)    (user)    (filter)    (bin)         (exp)     (object)
    regex = '([S_][F_][D_][B_])_([0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-?[0-9]{2}-[0-9]{2})_J([0-9]+)_+([A-Z]+)_([_A-Z0-9]+)_([0-9]x[0-9])_([0-9]+)s_(.*)'
    m = re.match(regex, filename)
    if m:
        flags, date, task_id, user, filter, bin, exp, object = m.groups()
        bin = int(bin[0])

        solved = (flags[0] == 'S')
        calibrated = (flags[1:4] == "FDB")

        # print(f"Parsed flags={flags} date={date} task_id={task_id} user={user} filter={filter} bin={bin} exp={exp} object={object}")
        return {
            "flags": flags,
            "date": date,
            "task_id": int(task_id),
            "user": user,
            "filter": filter,
            "binning": int(bin),
            "exposure": int(exp),
            "object": object,
            "imagename": fname,
            "solve": solved,
            "solved": solved,
            "calibrate": calibrated,
            "calibrated": calibrated
        }
    else:
        print(f"WARNING:Failed to parse filename {fname}")
        return {}
