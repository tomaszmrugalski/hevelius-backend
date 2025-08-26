"""
Commands for asteroid observation planning.
Handles downloading asteroid orbital data from Minor Planet Center,
caching it locally, and computing visibility for given observing conditions.
"""

import os
import sys
import requests
import sqlite3
from datetime import datetime, timedelta
from typing import List, Tuple, Optional
import math
import re

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.coordinates import solar_system_ephemeris
import astropy.coordinates as coord
import numpy as np

from hevelius import config


class AsteroidData:
    """Manages asteroid orbital data and computations."""

    def __init__(self, cache_dir: str = None):
        """Initialize with optional cache directory."""
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".hevelius", "asteroid_cache")

        self.cache_dir = cache_dir
        self.db_path = os.path.join(cache_dir, "asteroids.db")

        # Create cache directory if it doesn't exist
        os.makedirs(cache_dir, exist_ok=True)

        # Initialize local database
        self._init_db()

    def _init_db(self):
        """Initialize the local SQLite database for caching asteroid data."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Create asteroids table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS asteroids (
                number INTEGER PRIMARY KEY,
                name TEXT,
                epoch REAL,
                mean_anomaly REAL,
                perihelion_arg REAL,
                ascending_node REAL,
                inclination REAL,
                eccentricity REAL,
                mean_motion REAL,
                semimajor_axis REAL,
                absolute_magnitude REAL,
                slope_parameter REAL,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')

        # Create cache metadata table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS cache_metadata (
                key TEXT PRIMARY KEY,
                value TEXT,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')

        conn.commit()
        conn.close()

    def download_mpc_data(self, force_update: bool = False) -> bool:
        """
        Download asteroid orbital elements from Minor Planet Center.
        Returns True if successful, False otherwise.
        """
        # Check if we need to update (data older than 7 days)
        if not force_update and self._is_cache_fresh():
            print("Asteroid data is up to date (less than 7 days old)")
            return True

        print("Downloading asteroid orbital data from Minor Planet Center...")

        # MPC format: https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT
        mpc_url = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT"

        try:
            response = requests.get(mpc_url, timeout=300)  # 5 minute timeout
            response.raise_for_status()

            # Parse and store the data
            self._parse_mpc_data(response.text)

            # Update cache timestamp
            self._update_cache_timestamp()

            print(f"Successfully downloaded and cached asteroid data")
            return True

        except requests.RequestException as e:
            print(f"Error downloading asteroid data: {e}")
            return False
        except Exception as e:
            print(f"Error processing asteroid data: {e}")
            return False

    def _is_cache_fresh(self, max_age_days: int = 7) -> bool:
        """Check if cached data is fresh enough."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("SELECT value FROM cache_metadata WHERE key = 'last_update'")
        result = cursor.fetchone()
        conn.close()

        if not result:
            return False

        last_update = datetime.fromisoformat(result[0])
        age = datetime.now() - last_update

        return age.days < max_age_days

    def _update_cache_timestamp(self):
        """Update the cache timestamp."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute(
            "INSERT OR REPLACE INTO cache_metadata (key, value) VALUES (?, ?)",
            ('last_update', datetime.now().isoformat())
        )

        conn.commit()
        conn.close()

    def _parse_mpc_data(self, mpc_text: str):
        """Parse MPC orbital elements file and store in database."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Clear existing data
        cursor.execute("DELETE FROM asteroids")

        lines = mpc_text.strip().split('\n')
        count = 0

        for line in lines:
            # Skip header lines, empty lines, and separator lines
            if (len(line) < 200 or line.startswith('#') or not line.strip() or
                line.startswith('Des') or line.startswith('---') or
                line.startswith('MINOR PLANET') or line.startswith('This file') or
                line.strip() == ''):
                continue

            try:
                # Parse MPC format based on the actual format seen in the data
                # The format appears to be different from the documentation

                # Extract number from the beginning of the line
                number_str = line[0:5].strip()
                if not number_str or not number_str.isdigit():
                    continue

                number = int(number_str)

                # Extract name from the end of the line (after the last parentheses)
                name_match = re.search(r'\((\d+)\)\s+(.+?)(?:\s+\d{8})?$', line)
                if name_match:
                    name = name_match.group(2).strip()
                else:
                    name = ""

                # Parse the orbital elements from fixed positions
                # Based on the actual format observed
                parts = line.split()
                if len(parts) < 15:
                    continue

                # H magnitude (absolute magnitude)
                absolute_magnitude = float(parts[1]) if len(parts) > 1 else 99.0

                # G slope parameter
                slope_parameter = float(parts[2]) if len(parts) > 2 else 0.15

                # Epoch (packed format like K2555)
                epoch_str = parts[3] if len(parts) > 3 else ""
                # Convert packed epoch - this is a simplified conversion
                if epoch_str.startswith('K'):
                    # K2555 means 2025.55 (approximately)
                    year_part = epoch_str[1:]
                    if len(year_part) == 4:
                        year = 2000 + int(year_part[:2])
                        fraction = int(year_part[2:]) / 100.0
                        epoch = year + fraction
                    else:
                        epoch = 2025.0  # Default
                else:
                    epoch = 2025.0  # Default

                # Mean anomaly
                mean_anomaly = float(parts[4]) if len(parts) > 4 else 0.0

                # Argument of perihelion
                perihelion_arg = float(parts[5]) if len(parts) > 5 else 0.0

                # Longitude of ascending node
                ascending_node = float(parts[6]) if len(parts) > 6 else 0.0

                # Inclination
                inclination = float(parts[7]) if len(parts) > 7 else 0.0

                # Eccentricity
                eccentricity = float(parts[8]) if len(parts) > 8 else 0.0

                # Mean motion (daily motion)
                mean_motion = float(parts[9]) if len(parts) > 9 else 0.0

                # Semimajor axis
                semimajor_axis = float(parts[10]) if len(parts) > 10 else 0.0

                cursor.execute('''
                    INSERT INTO asteroids (
                        number, name, epoch, mean_anomaly, perihelion_arg,
                        ascending_node, inclination, eccentricity, mean_motion,
                        semimajor_axis, absolute_magnitude, slope_parameter
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    number, name, epoch, mean_anomaly, perihelion_arg,
                    ascending_node, inclination, eccentricity, mean_motion,
                    semimajor_axis, absolute_magnitude, slope_parameter
                ))

                count += 1

                if count % 10000 == 0:
                    print(f"Processed {count} asteroids...")

            except (ValueError, IndexError):
                # Skip malformed lines
                continue

        conn.commit()
        conn.close()

        print(f"Successfully processed {count} asteroids")

    def compute_visibility(self, location: EarthLocation, obs_date: str,
                          mag_min: float = 8.0, mag_max: float = 16.0,
                          alt_min: float = 20.0, constraint: str = None,
                          order_by: str = None) -> List[dict]:
        """
        Compute asteroid visibility for given location and date.

        Args:
            location: Observer location
            obs_date: Observation date (YYYY-MM-DD)
            mag_min: Minimum magnitude filter
            mag_max: Maximum magnitude filter
            alt_min: Minimum altitude filter (degrees)
            constraint: Additional SQL constraint (e.g., "number < 3000")
            order_by: Order by field (default: "absolute_magnitude")

        Returns:
            List of visible asteroids with their properties
        """

        # Convert date to astropy Time
        obs_time = Time(obs_date + ' 00:00:00')

        # Get asteroids from database
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Build SQL query with optional constraint and order_by
        base_query = """
            SELECT number, name, epoch, mean_anomaly, perihelion_arg,
                   ascending_node, inclination, eccentricity, mean_motion,
                   semimajor_axis, absolute_magnitude, slope_parameter
            FROM asteroids
            WHERE absolute_magnitude BETWEEN ? AND ?
        """

        params = [mag_min - 5, mag_max + 5]  # Add margin for apparent magnitude calculation

        # Add constraint if specified
        if constraint:
            base_query += f" AND ({constraint})"

        # Add ordering
        if order_by:
            # Validate order_by field to prevent SQL injection
            valid_fields = ['number', 'name', 'absolute_magnitude', 'semimajor_axis',
                           'eccentricity', 'inclination', 'mean_motion']
            if order_by in valid_fields:
                base_query += f" ORDER BY {order_by}"
            else:
                base_query += " ORDER BY absolute_magnitude"  # Default fallback
        else:
            base_query += " ORDER BY absolute_magnitude"  # Default ordering

        cursor.execute(base_query, params)
        asteroids = cursor.fetchall()
        conn.close()

        visible_asteroids = []

        print(f"Computing visibility for {len(asteroids)} asteroids...")

        # Define night time span (sunset to sunrise)
        night_times = self._get_night_times(location, obs_time)

        for i, asteroid in enumerate(asteroids):
            if i % 1000 == 0:
                print(f"Processed {i}/{len(asteroids)} asteroids...")

            try:
                # Compute position and magnitude at midnight
                midnight = obs_time + 0.5 * u.day  # Approximate midnight

                ra, dec, distance, magnitude = self._compute_asteroid_position(
                    asteroid, midnight
                )

                if magnitude < mag_min or magnitude > mag_max:
                    continue

                # Create SkyCoord for the asteroid
                asteroid_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

                # Compute altitude/azimuth for the night
                max_alt, max_alt_time = self._compute_max_altitude(
                    asteroid_coord, location, night_times
                )

                if max_alt < alt_min:
                    continue

                visible_asteroids.append({
                    'number': asteroid[0],
                    'name': asteroid[1],
                    'ra': ra,
                    'dec': dec,
                    'magnitude': magnitude,
                    'distance_au': distance,
                    'max_altitude': max_alt,
                    'max_alt_time': max_alt_time.iso if max_alt_time else None
                })

            except Exception:
                # Skip asteroids with computation errors
                continue

        print(f"Found {len(visible_asteroids)} visible asteroids")
        return visible_asteroids

    def _get_night_times(self, location: EarthLocation, obs_date: Time) -> List[Time]:
        """Get list of times throughout the night for calculations."""

        # Create time array from sunset to sunrise
        times = []
        for hour in range(24):
            time = obs_date + hour * u.hour
            times.append(time)

        return times

    def _compute_asteroid_position(self, asteroid_data: tuple, obs_time: Time) -> Tuple[float, float, float, float]:
        """
        Compute asteroid RA, Dec, distance, and apparent magnitude.

        This is a simplified calculation. For production use, consider
        using more sophisticated ephemeris calculations.
        """

        # Extract orbital elements
        (number, name, epoch, mean_anomaly, perihelion_arg,
         ascending_node, inclination, eccentricity, mean_motion,
         semimajor_axis, absolute_magnitude, slope_parameter) = asteroid_data

        # Convert epoch to JD
        epoch_jd = epoch

        # Time since epoch in days
        dt = obs_time.jd - epoch_jd

        # Mean anomaly at observation time
        M = mean_anomaly + mean_motion * dt
        M = M % 360.0

        # Solve Kepler's equation (simplified)
        E = M + (180.0/math.pi) * eccentricity * math.sin(math.radians(M))

        # True anomaly
        nu = 2 * math.atan(math.sqrt((1+eccentricity)/(1-eccentricity)) *
                          math.tan(math.radians(E/2)))
        nu = math.degrees(nu)

        # Distance from Sun
        r = semimajor_axis * (1 - eccentricity**2) / (1 + eccentricity * math.cos(math.radians(nu)))

        # Position in orbital plane
        x_orb = r * math.cos(math.radians(nu))
        y_orb = r * math.sin(math.radians(nu))

        # Convert to heliocentric ecliptic coordinates
        # This is simplified - proper implementation would use rotation matrices

        # For now, use approximate RA/Dec (this should be improved)
        # This is a placeholder calculation
        ra = (ascending_node + perihelion_arg + nu) % 360.0
        dec = inclination * math.sin(math.radians(nu))

        # Distance from Earth (approximation)
        distance = abs(r - 1.0)  # Rough approximation

        # Apparent magnitude calculation
        # V = H + 5*log10(r*delta) - 2.5*log10((1-G)*phi1 + G*phi2)
        # Simplified version:
        magnitude = absolute_magnitude + 5 * math.log10(r * distance)

        return ra, dec, distance, magnitude

    def _compute_max_altitude(self, coord: SkyCoord, location: EarthLocation,
                            times: List[Time]) -> Tuple[float, Optional[Time]]:
        """Compute maximum altitude during the night."""

        max_alt = -90.0
        max_alt_time = None

        for time in times:
            altaz = coord.transform_to(AltAz(obstime=time, location=location))
            alt = altaz.alt.degree

            if alt > max_alt:
                max_alt = alt
                max_alt_time = time

        return max_alt, max_alt_time


def asteroid_update(args):
    """Update asteroid orbital data from Minor Planet Center."""

    asteroid_data = AsteroidData()

    force = getattr(args, 'force', False)
    success = asteroid_data.download_mpc_data(force_update=force)

    if success:
        print("Asteroid data update completed successfully")
        return 0
    else:
        print("Failed to update asteroid data")
        return 1


def asteroid_visible(args):
    """List visible asteroids for given observing conditions."""

    # Parse location
    if hasattr(args, 'lat') and hasattr(args, 'lon'):
        location = EarthLocation(lat=args.lat*u.deg, lon=args.lon*u.deg,
                               height=args.elevation*u.m)
    else:
        # Default location (can be made configurable)
        location = EarthLocation(lat=50.0*u.deg, lon=20.0*u.deg, height=200*u.m)
        print(f"Using default location: lat=50°, lon=20°, elevation=200m")

    asteroid_data = AsteroidData()

    # Check if we have data
    if not os.path.exists(asteroid_data.db_path):
        print("No asteroid data found. Please run 'hevelius asteroid update' first.")
        return 1

    visible = asteroid_data.compute_visibility(
        location=location,
        obs_date=args.date,
        mag_min=args.mag_min,
        mag_max=args.mag_max,
        alt_min=args.alt_min,
        constraint=getattr(args, 'constraint', None),
        order_by=getattr(args, 'order_by', None)
    )

    # Print results
    if args.format == 'csv':
        print("Number,Name,RA,Dec,Magnitude,Distance_AU,Max_Altitude,Max_Alt_Time")
        for ast in visible:
            print(f"{ast['number']},{ast['name']},{ast['ra']:.4f},{ast['dec']:.4f},"
                  f"{ast['magnitude']:.2f},{ast['distance_au']:.3f},"
                  f"{ast['max_altitude']:.1f},{ast['max_alt_time']}")
    else:
        print(f"\nVisible asteroids for {args.date}:")
        print(f"Location: {location}")
        print(f"Magnitude range: {args.mag_min} - {args.mag_max}")
        print(f"Minimum altitude: {args.alt_min}°")
        print("-" * 80)
        print(f"{'Num':<8} {'Name':<20} {'RA':<10} {'Dec':<10} {'Mag':<6} {'Alt':<6} {'Time (UTC)'}")
        print("-" * 80)

        for ast in visible:
            print(f"{ast['number']:<8} {ast['name']:<20} "
                  f"{ast['ra']:<10.4f} {ast['dec']:<10.4f} "
                  f"{ast['magnitude']:<6.2f} {ast['max_altitude']:<6.1f} "
                  f"{ast['max_alt_time'] or 'N/A'}")

    print(f"\nTotal visible asteroids: {len(visible)}")
    return 0


def asteroid_info(args):
    """Show information about a specific asteroid."""

    asteroid_data = AsteroidData()

    if not os.path.exists(asteroid_data.db_path):
        print("No asteroid data found. Please run 'hevelius asteroid update' first.")
        return 1

    conn = sqlite3.connect(asteroid_data.db_path)
    cursor = conn.cursor()

    # Search by number or name
    if args.asteroid.isdigit():
        cursor.execute("SELECT * FROM asteroids WHERE number = ?", (int(args.asteroid),))
    else:
        cursor.execute("SELECT * FROM asteroids WHERE name LIKE ?", (f"%{args.asteroid}%",))

    results = cursor.fetchall()
    conn.close()

    if not results:
        print(f"Asteroid '{args.asteroid}' not found")
        return 1

    for asteroid in results[:10]:  # Limit to first 10 results
        print(f"\nAsteroid {asteroid[0]} ({asteroid[1]})")
        print(f"  Absolute magnitude: {asteroid[10]:.2f}")
        print(f"  Semimajor axis: {asteroid[9]:.3f} AU")
        print(f"  Eccentricity: {asteroid[7]:.3f}")
        print(f"  Inclination: {asteroid[6]:.2f}°")
        print(f"  Epoch: {asteroid[2]}")

    if len(results) > 10:
        print(f"\n... and {len(results) - 10} more matches")

    return 0