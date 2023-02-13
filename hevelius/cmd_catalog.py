from hevelius import db
from hevelius.utils import parse_dec, parse_ra, format_dec, format_ra
from hevelius import config
from argparse import ArgumentTypeError
import sys


def format_get(format: str) -> str:
    """
    Ensures that the specified output format allowed. Allowed values are:
    - none - don't print frames list at all
    - filenames - print filenames only (useful for scripting)
    - csv - print all found frames in CSV format
    - brief - print a space separated list on the screen
    - full - print everything on the screen
    - pixinsight - print a list of frames in PixInsight format (useful for
                    importing into SubframeSelector)
    """
    if format not in ["none", "filenames", "csv", "brief", "full", "pixinsight"]:
        raise ArgumentTypeError(f"unsupported format: {format}, allowed are: none, filenames, csv, brief, full, pixinsight")

    return format


def catalog(args):
    """
    Searches the objects catalog for objects close to specified coordinates.
    """

    conn = db.connect()

    ra = 0
    decl = 0
    if len(args.object):
        print(f"Looking for object {args.object}")
        obj = db.catalog_get(conn, args.object)
        if obj != []:
            print(obj)
            ra = obj[0][3]
            decl = obj[0][4]
            print(f"Found object {args.object} in a catalog, using its coords: RA {ra} DEC {decl}")
    else:
        ra = parse_ra(args.ra)
        decl = parse_dec(args.decl)

    radius = float(args.proximity)
    format = format_get(args.format)

    if format != "csv":
        print(f"Looking for objects close to RA {format_ra(ra)} DEC {format_dec(decl)}, within a radius of {radius} deg")

        objects = db.catalog_radius_get(conn, ra, decl, radius)

        print(f"Found {len(objects)} object(s) that match criteria: distance from RA {format_ra(ra)} DEC {format_dec(decl)} no larger than {radius}")
        for object in objects:
            print(object)

    frames = db.tasks_radius_get(conn, ra, decl, radius)

    if format != "csv":
        print(f"Found {len(frames)} frame(s) that match criteria: distance from RA {format_ra(ra)} DEC {format_dec(decl)} no larger than {radius}")
    conn.close()

    if format == "none":
        return

    if format == "csv":
        print("# task_id, object, filename, fwhm, ra, decl, comment")
    # Print a space separated list of task IDs
    for frame in frames:
        if format == "filenames":
            print(frame[2])
        elif format == "brief":
            sys.stdout.write(f"{frame[0]} ")
        elif format == "full":
            print(f"Task {frame[0]}: RA {frame[4]} DEC {frame[5]}, object: {frame[1]}, file: {frame[2]}, fwhm: {frame[3]}")
        elif format == "csv":
            print(f"{frame[0]},{frame[1]},{frame[2]},{frame[3]},{frame[4]},{frame[5]},{frame[6]}")
        elif format == "pixinsight":
            print(f'   [true, "{config.REPO_PATH}/{frame[2]}", "", ""],')
    print("")
