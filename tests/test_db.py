from hevelius import db
import unittest
from tests.dbtest import use_repository


class DbTest(unittest.TestCase):

    @use_repository
    def test_db_version_get(self, config):
        """Test that the database schema version is correct."""

        conn = db.connect(config)
        version = db.version_get(conn)
        conn.close()

        self.assertEqual(version, 11)

    @use_repository
    def test_sensor(self, config):
        """Test that the sensors information can be retrieved."""
        conn = db.connect(config)

        cases = [
            {
                "name": "FLI",
                "exp_sensor_id": 2,
                "exp_name": "FLI Proline 16803",
                "exp_resx": 4096,
                "exp_resy": 4096,
                "exp_pixel_x": 9.0,
                "exp_pixel_y": 9.0,
                "exp_bits": 16,
                "exp_width": 36.8,
                "exp_height": 36.8
            }
        ]

        for case in cases:
            sensor_id, name, resx, resy, pixel_x, pixel_y, bits, width, height = db.sensor_get_by_name(conn, case["name"])
            self.assertAlmostEqual(sensor_id, case['exp_sensor_id'])
            self.assertAlmostEqual(name, case['exp_name'])
            self.assertAlmostEqual(resx, case['exp_resx'])
            self.assertAlmostEqual(resy, case['exp_resy'])
            self.assertAlmostEqual(pixel_x, case['exp_pixel_x'])
            self.assertAlmostEqual(pixel_y, case['exp_pixel_y'])
            self.assertAlmostEqual(bits, case['exp_bits'])
            self.assertAlmostEqual(width, case['exp_width'])
            self.assertAlmostEqual(height, case['exp_height'])

        conn.close()

    @use_repository
    def test_tasks_radius_get(self, config):
        """Test that tasks can be retrieved by filter."""

        conn = db.connect(config)
        tasks = db.tasks_radius_get(conn, 0, 25, 1, "", "task_id DESC")
        conn.close()

        self.assertGreaterEqual(len(tasks), 10)
        self.assertEqual(tasks[0][0], 87775)
        self.assertEqual(tasks[1][0], 70556)
        self.assertEqual(tasks[2][0], 69426)
        # Let's assume the remaining tasks have correct task_ids, too.

    @use_repository
    def test_tasks_radius_get_filter(self, config):
        """Test that tasks can be retrieved by filter."""

        cases = [
            {"filter": "AND scope_id=1", "exp_count": 7},
            {"filter": "AND scope_id=2", "exp_count": 4},
            {"filter": "AND scope_id=4", "exp_count": 0},
            {"filter": "AND sensor_id=1", "exp_count": 1},
            {"filter": "AND sensor_id=2", "exp_count": 10},
            {"filter": "AND he_resx=3326", "exp_count": 1},
            {"filter": "AND he_resy=2504", "exp_count": 1},
            {"filter": "AND scope_id=4", "exp_count": 0}
        ]

        conn = db.connect(config)

        for case in cases:
            tasks = db.tasks_radius_get(conn, 0, 0, 360, case['filter'], "")
            self.assertEqual(len(tasks), case['exp_count'])

        conn.close()
