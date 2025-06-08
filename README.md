# Hevelius Backend

![pylint](https://github.com/tomaszmrugalski/hevelius-backend/actions/workflows/pylint.yml/badge.svg)
![pytest](https://github.com/tomaszmrugalski/hevelius-backend/actions/workflows/testing.yml/badge.svg)
![CodeQL](https://github.com/tomaszmrugalski/hevelius-backend/actions/workflows/github-code-scanning/codeql/badge.svg)

Hevelius is an astronomical data processing and management system.

## Features

- Database management for astronomical observations
- File repository handling for FITS files
- Statistical analysis of observation data
- Catalog searches for astronomical objects
- Asteroid observation planning

## Installation

1. Clone the repository
2. Create and activate virtual environment:

   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

## Usage

The main command-line interface is provided through `bin/hevelius`.

### Basic Commands

- `hevelius version` - Show Hevelius version
- `hevelius config` - Show current configuration
- `hevelius db version` - Show database schema version
- `hevelius db migrate` - Migrate database to latest schema
- `hevelius db backup` - Create database backup
- `hevelius db stats` - Show database statistics

### Repository Management

- `hevelius repo -f <file>` - Process single FITS file
- `hevelius repo -l <list>` - Process files from list
- `hevelius repo -d <dir>` - Process all FITS files in directory recursively

### Data Analysis

- `hevelius data distrib` - Show photo distribution
- `hevelius data groups` - Show frame groups
- `hevelius data catalog` - Find astronomical objects in catalog

### Asteroid Observation Planning

The asteroid planning functionality helps astronomers plan observations of asteroids by:

1. **Downloading and caching orbital data** from the Minor Planet Center
2. **Computing visibility** for specific locations and dates
3. **Filtering by magnitude and altitude** to find suitable targets
4. **Providing detailed information** about specific asteroids

#### Update Asteroid Data

Download the latest asteroid orbital elements from the Minor Planet Center:

```bash
hevelius asteroid update
```

Use `--force` to update even if data is recent:

```bash
hevelius asteroid update --force
```

#### Find Visible Asteroids

List asteroids visible for a specific date and location:

```bash
hevelius asteroid visible --date 2024-12-25 --lat 50.0 --lon 20.0
```

**Parameters:**

- `--date` (required): Observation date in YYYY-MM-DD format
- `--lat`: Observer latitude in degrees (default: 50.0)
- `--lon`: Observer longitude in degrees (default: 20.0)
- `--elevation`: Observer elevation in meters (default: 200.0)
- `--mag-min`: Minimum magnitude filter (default: 8.0)
- `--mag-max`: Maximum magnitude filter (default: 16.0)
- `--alt-min`: Minimum altitude filter in degrees (default: 20.0)
- `--format`: Output format - `table` or `csv` (default: table)

**Example with filters:**

```bash
hevelius asteroid visible --date 2024-12-25 --lat 52.5 --lon 13.4 \
    --mag-min 10.0 --mag-max 14.0 --alt-min 30.0 --format csv
```

This will find asteroids between magnitude 10-14 that reach at least 30° altitude for Berlin on Christmas Day 2024.

#### Get Asteroid Information

Search for information about a specific asteroid by number or name:

```bash
hevelius asteroid info 1
hevelius asteroid info Ceres
hevelius asteroid info Vesta
```

This is a backend interface for Hevelius, an astronomy processing software and
observatory management system. It's in the early stages of development, but some
of the features are usable already.

## Current capabilities (command-line)

Status as of June 2025:

- **Scan FITS repository on disk**: Hevelius is able to scan a local disk for FITS files, extract some data from found files
  (from filenames and FITS header) and put this information into PostgreSQL DB. Then the DB is used to report various
  characteristics.
- **Full sky histogram**: Generate full sky distribution of found frames with 1 degree resolution. The data is presented as
  interactive RA/DEC chart.
- **Points of interest**: Generate a list of the most commonly photographed coordinates in the sky.
- **Objects and frames search**: Ability to find catalog objects and frames based on specified RA/DEC coordinates and radius.
- **4 Catalogs**: Provides 4 catalogs (NGC, IC, Messier, and Caldwell) in a DB format and a basic interface to query it.
- **PixInsight integration**: This is in the very early stages. The idea is that Hevelius will be able to offload certain
  tasks to Pix or at least export/import data in a format that's compatible with PixInsight.
- **Command line interface**: Currently Hevelius has a command line interface written in `python`. A Rest API and gui front-end
  is planned, but currently not a priority.
- **Ability to search based on distance**. Implemented proper Haversine formula.
- **Database management**: Schema versioning and upgrades, backup, etc.
- **Configuration**: Config file support and some limited environment variables.
- **Asteroid visibility**: Retrieves orbital data from MPC (Minor Planets Center) and lists asteroid that are locally visible


### Asteroid visibility

The visibility command provides:

- **Asteroid number and name**
- **Right Ascension and Declination** (J2000.0)
- **Apparent magnitude** at observation date
- **Maximum altitude** during the night
- **Time of maximum altitude** (UTC)

Example output:
```
Visible asteroids for 2024-12-25:
Location: <EarthLocation (1234567. m, 5678901. m, 2345678. m) in geocentric coords>
Magnitude range: 10.0 - 14.0
Minimum altitude: 30.0°
--------------------------------------------------------------------------------
Num      Name                 RA         Dec        Mag    Alt    Time (UTC)
--------------------------------------------------------------------------------
4        Vesta               123.4567   -12.3456   12.34  45.6   2024-12-25T23:45:00.000
7        Iris                234.5678   +23.4567   13.21  38.9   2024-12-25T22:15:00.000
8        Flora               345.6789   -5.6789    11.87  52.3   2024-12-26T01:30:00.000

Total visible asteroids: 3
```

## Current capabilities (REST API)

- Log in users
- List tasks (with pagination, filtering, and sorting)
- Add new observation task
- Edit existing observation task
- List objects from catalogs (Messier, Caldwell, NGC, IC are supported)
- Search objects in catalogs
- Display heat map of the sky (sky map colored with number of photos taken in each square degree)

## Documentation

- [Installation](doc/install.md) - You probably want to start here.
- [Commands reference](doc/commands.md) - Available commands are (or soon will) be documented here.
- [Catalogs](doc/catalogs.md) - Hevelius comes with several astronomical catalogs.
- [Database details](doc/db.md) - The most useful section is probably the paragraph about DB initalization.

# Dependencies

Key dependencies include:
- `astropy` - Astronomical calculations and coordinate transformations
- `psycopg2-binary` - PostgreSQL database interface
- `pandas` - Data analysis and manipulation
- `numpy` - Numerical computations
- `requests` - HTTP requests for downloading asteroid data
- `flask` - Web framework for API endpoints

## Developer's corner

- [Developer's guide](doc/devel.md)
- [Security info](SECURITY.md)
- [License](LICENSE)
