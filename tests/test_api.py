# This is a workaround to suppress the specific marshmallow warning
# that happens when running the tests. This warning is not relevant
# and it's caused by one of our dependencies.

import unittest
import os
import warnings
import json
from datetime import datetime
from flask_jwt_extended import create_access_token
from tests.dbtest import use_repository
from marshmallow import warnings as marshmallow_warnings
from hevelius import db

# Suppress the specific marshmallow warning
warnings.filterwarnings("ignore", category=marshmallow_warnings.RemovedInMarshmallow4Warning)

from heveliusbackend.app import app  # noqa: E402


class TestTaskAdd(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        self.app.testing = True

        # Create a test JWT token
        with app.app_context():
            self.test_token = create_access_token(
                identity=1,  # user_id=1
                additional_claims={
                    'permissions': 1,
                    'username': 'test_user'
                }
            )
            self.headers = {
                'Authorization': f'Bearer {self.test_token}',
                'Content-Type': 'application/json'
            }

    @use_repository
    def test_task_add_success(self, config):
        """Test successful task addition"""
        test_task = {
            "user_id": 1,  # This should match the token's identity
            "scope_id": 1,  # Required telescope ID
            "object": "M31",
            "ra": 0.712,  # ~00h 42m for M31
            "decl": 41.27,  # ~41° 16' for M31
            "exposure": 300.0,
            "filter": "L",
            "binning": 1,
            "guiding": True,
            "dither": False,
            "solve": True,
            "calibrate": True
        }

        os.environ['HEVELIUS_DB_NAME'] = config['database']

        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertTrue(data['status'])
        self.assertIsInstance(data['task_id'], int)
        self.assertIn('Task', data['msg'])
        self.assertIn('created successfully', data['msg'])

        # Verify the task was actually added to the database
        task_id = data['task_id']
        cnx = db.connect()
        query = """SELECT task_id, user_id, scope_id, object, ra, decl, exposure,
                         filter, binning, guiding, dither, solve, calibrate, state
                  FROM tasks WHERE task_id = %s"""
        result = db.run_query(cnx, query, (task_id,))
        cnx.close()

        self.assertIsNotNone(result, "Task not found in database")
        self.assertEqual(len(result), 1, "Expected exactly one task")

        task = result[0]
        # Verify all fields match what we sent
        self.assertEqual(task[0], task_id)  # task_id
        self.assertEqual(task[1], test_task['user_id'])  # user_id
        self.assertEqual(task[2], test_task['scope_id'])  # scope_id
        self.assertEqual(task[3], test_task['object'])  # object
        self.assertEqual(float(task[4]), test_task['ra'])  # ra
        self.assertEqual(float(task[5]), test_task['decl'])  # decl
        self.assertEqual(float(task[6]), test_task['exposure'])  # exposure
        self.assertEqual(task[7], test_task['filter'])  # filter
        self.assertEqual(task[8], test_task['binning'])  # binning
        self.assertEqual(bool(task[9]), test_task['guiding'])  # guiding
        self.assertEqual(bool(task[10]), test_task['dither'])  # dither
        self.assertEqual(bool(task[11]), test_task['solve'])  # solve
        self.assertEqual(bool(task[12]), test_task['calibrate'])  # calibrate
        self.assertEqual(task[13], 1)  # state should be 1 for new tasks

        os.environ.pop('HEVELIUS_DB_NAME')

    def test_task_add_missing_required(self):
        """Test task addition with missing required fields"""
        test_task = {
            "object": "M31",
            "exposure": 300.0
            # Missing required fields: user_id and scope_id
        }

        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 422)  # Unprocessable Entity
        self.assertIn('errors', data)

    def test_task_add_invalid_ra(self):
        """Test task addition with invalid RA value"""
        test_task = {
            "user_id": 1,
            "scope_id": 1,
            "object": "M31",
            "ra": 25.0,  # Invalid: RA must be 0-24
            "decl": 41.27,
            "exposure": 300.0
        }

        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 422)  # Unprocessable Entity
        self.assertIn('errors', data)

    def test_task_add_invalid_dec(self):
        """Test task addition with invalid declination value"""
        test_task = {
            "user_id": 1,
            "scope_id": 1,
            "object": "M31",
            "ra": 0.712,
            "decl": 91.0,  # Invalid: Dec must be -90 to +90
            "exposure": 300.0
        }

        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 422)  # Unprocessable Entity
        self.assertIn('errors', data)

    def test_task_add_unauthorized(self):
        """Test task addition with mismatched user_id"""
        test_task = {
            "user_id": 2,  # Different from token's identity (1)
            "scope_id": 1,
            "object": "M31",
            "ra": 0.712,
            "decl": 41.27,
            "exposure": 300.0
        }

        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)  # API returns 200 with error message
        self.assertFalse(data['status'])
        self.assertEqual(data['msg'], 'Unauthorized: token user_id does not match request user_id')


class TestVersion(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        self.app.testing = True

    def test_version_endpoint(self):
        """Test version endpoint returns correct version"""
        from hevelius.version import VERSION

        response = self.app.get('/api/version')
        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertIn('version', data)
        self.assertEqual(data['version'], VERSION)


class TestTaskGet(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        self.app.testing = True

        # Create a test JWT token
        with app.app_context():
            self.test_token = create_access_token(
                identity=1,  # user_id=1
                additional_claims={
                    'permissions': 1,
                    'username': 'test_user'
                }
            )
            self.headers = {
                'Authorization': f'Bearer {self.test_token}',
                'Content-Type': 'application/json'
            }

    @use_repository
    def test_task_get_success(self, config):
        """Test successful task retrieval"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # First create a task
        test_task = {
            "user_id": 1,
            "scope_id": 1,  # Required telescope ID
            "object": "M31",
            "ra": 0.712,
            "decl": 41.27,
            "exposure": 300.0,
            "filter": "L",
            "binning": 1,
            "guiding": True,
            "dither": False,
            "solve": True,
            "calibrate": True
        }

        # Add the task
        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)

        data = json.loads(response.data)
        task_id = data['task_id']

        # Now try to get the task
        response = self.app.get(f'/api/task-get?task_id={task_id}',
                                headers=self.headers)

        os.environ.pop('HEVELIUS_DB_NAME')

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertTrue(data['status'])
        self.assertIsNotNone(data['task'])
        self.assertEqual(data['task']['object'], 'M31')
        self.assertEqual(data['task']['ra'], 0.712)
        self.assertEqual(data['task']['decl'], 41.27)
        self.assertTrue(data['task']['scope_id'])
        self.assertTrue(test_task['scope_id'])
        self.assertEqual(data['task']['scope_id'], test_task['scope_id'])

    @use_repository
    def test_task_get_not_found(self, config):
        """Test task retrieval with non-existent task ID"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        response = self.app.get('/api/task-get?task_id=999999',
                                headers=self.headers)

        os.environ.pop('HEVELIUS_DB_NAME')

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertFalse(data['status'])
        self.assertIn('not found', data['msg'].lower())
        self.assertIsNone(data['task'])

    def test_task_get_no_auth(self):
        """Test task retrieval without authentication"""
        response = self.app.get('/api/task-get?task_id=1')
        self.assertEqual(response.status_code, 401)  # Unauthorized

    def test_task_get_missing_id(self):
        """Test task retrieval without task_id parameter"""
        response = self.app.get('/api/task-get',
                                headers=self.headers)
        self.assertEqual(response.status_code, 422)  # Unprocessable Entity


class TestTaskUpdate(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        self.app.testing = True

        # Create a test JWT token
        with app.app_context():
            self.test_token = create_access_token(
                identity=1,  # user_id=1
                additional_claims={
                    'permissions': 1,
                    'username': 'test_user'
                }
            )
            self.headers = {
                'Authorization': f'Bearer {self.test_token}',
                'Content-Type': 'application/json'
            }

    @use_repository
    def test_task_update_success(self, config):
        """Test successful task update"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # First create a task
        test_task = {
            "user_id": 1,
            "scope_id": 1,  # Required telescope ID
            "object": "M31",
            "ra": 0.712,
            "decl": 41.27,
            "exposure": 300.0,
            "filter": "L",
            "binning": 1,
            "guiding": True,
            "dither": False,
            "solve": True,
            "calibrate": True
        }

        # Add the task
        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)
        data = json.loads(response.data)
        task_id = data['task_id']

        # Update the task
        update_data = {
            "task_id": task_id,
            "scope_id": 2,  # Test updating scope_id
            "object": "M33",
            "exposure": 600.0
        }

        response = self.app.post('/api/task-update',
                                 data=json.dumps(update_data),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertTrue(data['status'])
        self.assertIn('updated successfully', data['msg'])

        # Verify the update
        response = self.app.get(f'/api/task-get?task_id={task_id}',
                                headers=self.headers)
        data = json.loads(response.data)

        self.assertEqual(data['task']['object'], 'M33')
        self.assertEqual(data['task']['exposure'], 600.0)
        # Original fields should remain unchanged
        self.assertEqual(data['task']['ra'], 0.712)
        self.assertEqual(data['task']['decl'], 41.27)
        self.assertEqual(data['task']['scope_id'], 2)  # Verify scope_id was updated

        os.environ.pop('HEVELIUS_DB_NAME')

    @use_repository
    def test_task_update_not_found(self, config):
        """Test updating non-existent task"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        update_data = {
            "task_id": 999999,
            "object": "M33"
        }

        response = self.app.post('/api/task-update',
                                 data=json.dumps(update_data),
                                 headers=self.headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertFalse(data['status'])
        self.assertIn('not found', data['msg'].lower())

        os.environ.pop('HEVELIUS_DB_NAME')

    def test_task_update_no_auth(self):
        """Test task update without authentication"""
        update_data = {
            "task_id": 1,
            "object": "M33"
        }

        # Send request without headers (no authentication)
        response = self.app.post('/api/task-update',
                                 data=json.dumps(update_data))
        self.assertEqual(response.status_code, 401)  # Unauthorized

    def test_task_update_missing_task_id(self):
        """Test task update without task_id"""
        update_data = {
            "object": "M33"
        }

        response = self.app.post('/api/task-update',
                                 data=json.dumps(update_data),
                                 headers=self.headers)
        self.assertEqual(response.status_code, 422)  # Unprocessable Entity

    @use_repository
    def test_task_update_unauthorized_user(self, config):
        """Test updating task owned by different user"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # First create a task as user 1
        test_task = {
            "user_id": 1,
            "scope_id": 1,
            "object": "M31",
            "ra": 0.712,
            "decl": 41.27
        }

        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)
        task_id = json.loads(response.data)['task_id']

        # Create token for different user
        with app.app_context():
            other_token = create_access_token(
                identity=2,  # Different user_id
                additional_claims={
                    'permissions': 1,
                    'username': 'other_user'
                }
            )
            other_headers = {
                'Authorization': f'Bearer {other_token}',
                'Content-Type': 'application/json'
            }

        # Try to update the task as different user
        update_data = {
            "task_id": task_id,
            "object": "M33"
        }

        response = self.app.post('/api/task-update',
                                 data=json.dumps(update_data),
                                 headers=other_headers)

        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertFalse(data['status'])
        self.assertIn('unauthorized', data['msg'].lower())

        os.environ.pop('HEVELIUS_DB_NAME')


class TestNightPlan(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        app.testing = True  # Set testing flag on the actual app instance

        # Create a test JWT token
        with app.app_context():
            self.test_token = create_access_token(
                identity=1,  # user_id=1
                additional_claims={
                    'permissions': 1,
                    'username': 'test_user'
                }
            )
            self.headers = {
                'Authorization': f'Bearer {self.test_token}',
                'Content-Type': 'application/json'
            }

    def tearDown(self):
        """Clean up after each test"""
        app.testing = False  # Reset testing flag

    @use_repository
    def test_night_plan_success(self, config):
        """Test successful night plan retrieval"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # First create test users
        cnx = db.connect(config)
        db.run_query(cnx, """
            INSERT INTO users (user_id, login, pass, firstname, lastname, share, phone, email, permissions)
            VALUES
            (100, 'test_user', 'test_pass', 'Test', 'User', 1.0, '123456789', 'test@test.com', 1),
            (101, 'other_user', 'test_pass', 'Other', 'User', 1.0, '987654321', 'other@test.com', 1)
            RETURNING user_id""")  # Add RETURNING clause to make it return results
        cnx.close()

        # Now create test tasks
        test_tasks = [
            {
                "user_id": 100,
                "scope_id": 1,
                "object": "M31",
                "ra": 0.712,
                "decl": 41.27,
                "exposure": 300.0,
                "state": 1  # New task
            },
            {
                "user_id": 100,
                "scope_id": 1,
                "object": "M33",
                "ra": 1.5,
                "decl": 30.0,
                "exposure": 300.0,
                "state": 1  # Should NOT be included (Template task)
            },
            {
                "user_id": 101,  # Different user
                "scope_id": 1,
                "object": "M51",
                "ra": 13.5,
                "decl": 47.0,
                "exposure": 300.0,
                "state": 1  # Should be included when not filtering by user
            },
            {
                "user_id": 100,
                "scope_id": 2,  # Different scope
                "object": "M45",
                "ra": 3.75,
                "decl": 24.1,
                "exposure": 300.0,
                "state": 1  # Should not be included due to scope_id
            }
        ]

        # Add test tasks
        task_ids = []
        for task in test_tasks:
            response = self.app.post('/api/task-add',
                                     data=json.dumps(task),
                                     headers=self.headers)
            data = json.loads(response.data)
            self.assertTrue(data['status'])
            task_ids.append(data['task_id'])

        # Test night plan without user_id filter
        response = self.app.get('/api/night-plan?scope_id=1',
                                headers=self.headers)

        data = json.loads(response.data)
        self.assertEqual(response.status_code, 200)
        self.assertIn('tasks', data)
        tasks = data['tasks']

        # Should find 3 tasks (all tasks for scope_id=1)
        self.assertEqual(len(tasks), 3)
        task_objects = [task['object'] for task in tasks]
        self.assertIn('M31', task_objects)
        self.assertIn('M51', task_objects)

        # Test night plan with user_id filter
        response = self.app.get('/api/night-plan?scope_id=1&user_id=100',
                                headers=self.headers)

        data = json.loads(response.data)
        self.assertEqual(response.status_code, 200)
        self.assertIn('tasks', data)
        tasks = data['tasks']

        # Should find 2 tasks (only user_id=100 tasks for scope_id=1)
        self.assertEqual(len(tasks), 2)
        task_objects = [task['object'] for task in tasks]
        self.assertIn('M31', task_objects)
        self.assertNotIn('M51', task_objects)

        os.environ.pop('HEVELIUS_DB_NAME')

    def test_night_plan_no_auth(self):
        """Test night plan retrieval without authentication"""
        response = self.app.get('/api/night-plan?scope_id=1')
        self.assertEqual(response.status_code, 401)  # Unauthorized

    def test_night_plan_missing_scope(self):
        """Test night plan retrieval without scope_id"""
        response = self.app.get('/api/night-plan',
                                headers=self.headers)
        self.assertEqual(response.status_code, 422)  # Unprocessable Entity

    @use_repository
    def test_night_plan_with_date_filters(self, config):
        """Test night plan with date parameter and skip_before/skip_after dates"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        test_date = datetime(2024, 3, 15).date()  # Specific test date
        before_date = datetime(2024, 3, 14).date()
        after_date = datetime(2024, 3, 16).date()

        # Create test tasks with different date constraints
        test_tasks = [
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "M31",
                "ra": 0.712,
                "decl": 41.27,
                "exposure": 300.0,
                "skip_before": datetime(2024, 3, 14).isoformat(),  # Before test date
                "skip_after": datetime(2024, 3, 16).isoformat(),   # After test date
                "state": 1  # Should be included
            },
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "M33",
                "ra": 1.5,
                "decl": 30.0,
                "exposure": 300.0,
                "skip_before": datetime(2024, 3, 16).isoformat(),  # After test date
                "state": 1  # Should not be included
            },
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "M51",
                "ra": 13.5,
                "decl": 47.0,
                "exposure": 300.0,
                "skip_after": datetime(2024, 3, 14).isoformat(),   # Before test date
                "state": 1  # Should not be included
            },
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "M42",
                "ra": 5.5,
                "decl": -5.4,
                "exposure": 300.0,
                # No date constraints
                "state": 1  # Should be included
            }
        ]

        # Add test tasks
        for task in test_tasks:
            response = self.app.post('/api/task-add',
                                     data=json.dumps(task),
                                     headers=self.headers)
            self.assertTrue(json.loads(response.data)['status'])

        # Test 1: Get night plan for specific date
        response = self.app.get(
            f'/api/night-plan?scope_id=1&user_id=1&date={test_date.isoformat()}',
            headers=self.headers
        )

        data = json.loads(response.data)
        self.assertEqual(response.status_code, 200)
        self.assertIn('tasks', data)
        tasks = data['tasks']

        # Should find 2 tasks (M31 which is within range and M42 which has no constraints)
        self.assertEqual(len(tasks), 2)
        task_objects = set(task['object'] for task in tasks)
        self.assertEqual(task_objects, {'M31', 'M42'})

        # Test 2: Get night plan for date before skip_before
        response = self.app.get(
            f'/api/night-plan?scope_id=1&user_id=1&date={before_date.isoformat()}',
            headers=self.headers
        )

        data = json.loads(response.data)
        tasks = data['tasks']
        task_objects = set(task['object'] for task in tasks)
        # Should only find M42 (no date constraints)
        self.assertEqual(task_objects, {'M42'})

        # Test 3: Get night plan for date after skip_after
        response = self.app.get(
            f'/api/night-plan?scope_id=1&user_id=1&date={after_date.isoformat()}',
            headers=self.headers
        )

        data = json.loads(response.data)
        tasks = data['tasks']
        task_objects = set(task['object'] for task in tasks)
        # Should only find M42 (no date constraints)
        self.assertEqual(task_objects, {'M42'})

        os.environ.pop('HEVELIUS_DB_NAME')

    @use_repository
    def test_night_plan_invalid_date_format(self, config):
        """Test night plan with invalid date format"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # Test with invalid date format
        response = self.app.get(
            '/api/night-plan?scope_id=1&user_id=1&date=invalid-date',
            headers=self.headers
        )

        self.assertEqual(response.status_code, 422)  # Unprocessable Entity

        os.environ.pop('HEVELIUS_DB_NAME')


class TestTasks(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        app.testing = True  # Set testing flag on the actual app instance

        # Create a test JWT token
        with app.app_context():
            self.test_token = create_access_token(
                identity=1,  # user_id=1
                additional_claims={
                    'permissions': 1,
                    'username': 'test_user'
                }
            )
            self.headers = {
                'Authorization': f'Bearer {self.test_token}',
                'Content-Type': 'application/json'
            }

    @use_repository(load_test_data="tests/test-data-basic.psql")
    def test_tasks_pagination(self, config):
        """Test tasks endpoint pagination"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # Create 150 test tasks
        test_tasks = []
        for i in range(150):
            task = {
                "user_id": 1,
                "scope_id": 1,
                "object": f"Test Object {i}",
                "ra": 0.712,
                "decl": 41.27,
                "exposure": 300.0,
                "state": 1
            }
            test_tasks.append(task)

        # Add all test tasks
        for task in test_tasks:
            response = self.app.post('/api/task-add',
                                     data=json.dumps(task),
                                     headers=self.headers)
            self.assertTrue(json.loads(response.data)['status'])

        # Test default pagination (page 1, 100 per page)
        response = self.app.get('/api/tasks', headers=self.headers)
        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(data['tasks']), 100)  # Default per_page
        self.assertEqual(data['page'], 1)
        self.assertEqual(data['per_page'], 100)
        self.assertEqual(data['total'], 150)
        self.assertEqual(data['pages'], 2)

        # Test second page
        response = self.app.get('/api/tasks?page=2', headers=self.headers)
        data = json.loads(response.data)

        self.assertEqual(len(data['tasks']), 50)  # Remaining tasks
        self.assertEqual(data['page'], 2)

        # Test custom per_page
        response = self.app.get('/api/tasks?per_page=50', headers=self.headers)
        data = json.loads(response.data)

        self.assertEqual(len(data['tasks']), 50)
        self.assertEqual(data['pages'], 3)

        os.environ.pop('HEVELIUS_DB_NAME')

    @use_repository(load_test_data="tests/test-data-basic.psql")
    def test_tasks_sorting(self, config):
        """Test tasks endpoint sorting"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # Create test tasks with different values
        test_tasks = [
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "A Object",
                "ra": 1.0,
                "decl": 0.0,
                "exposure": 100.0,
                "state": 1
            },
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "B Object",
                "ra": 2.0,
                "decl": 10.0,
                "exposure": 200.0,
                "state": 0
            }
        ]

        # Add test tasks
        for task in test_tasks:
            response = self.app.post('/api/task-add',
                                     data=json.dumps(task),
                                     headers=self.headers)
            self.assertTrue(json.loads(response.data)['status'])

        # Test sorting by different fields
        sort_tests = [
            ('object', 'asc', 'A Object'),
            ('object', 'desc', 'B Object'),
            ('ra', 'asc', 1.0),
            ('ra', 'desc', 2.0),
            ('exposure', 'asc', 100.0),
            ('state', 'asc', 0)
        ]

        # TODO: Sorting by state doesn't seem to work.

        for sort_by, sort_order, expected_first in sort_tests:
            response = self.app.get(
                f'/api/tasks?sort_by={sort_by}&sort_order={sort_order}',
                headers=self.headers
            )
            data = json.loads(response.data)

            self.assertEqual(response.status_code, 200)
            self.assertTrue(len(data['tasks']) > 0)

            # Check if sorting worked
            first_task = data['tasks'][0]
            self.assertEqual(first_task[sort_by], expected_first)

        os.environ.pop('HEVELIUS_DB_NAME')

    @use_repository(load_test_data="tests/test-data-basic.psql")
    def test_tasks_filtering(self, config):
        """Test tasks endpoint filtering"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # Create test tasks with different values
        test_tasks = [
            {
                "user_id": 1,
                "scope_id": 1,
                "object": "M31",
                "ra": 0.712,
                "decl": 41.27,
                "exposure": 300.0,
                "state": 1,
                "descr": "Test description 1"
            },
            {
                "user_id": 2,
                "scope_id": 2,
                "object": "M33",
                "ra": 1.5,
                "decl": 30.0,
                "exposure": 400.0,
                "state": 0,
                "descr": "Test description 2"
            }
        ]

        # Add test tasks
        for task in test_tasks:
            response = self.app.post('/api/task-add',
                                     data=json.dumps(task),
                                     headers=self.headers)
            self.assertTrue(json.loads(response.data)['status'])

        # Test various filters
        filter_tests = [
            ('user_id=1', 1),
            ('scope_id=2', 1),
            ('object=M31', 1),
            ('ra_min=1.0&ra_max=2.0', 1),
            ('decl_min=35.0&decl_max=45.0', 1),
            ('exposure=300.0', 1),
            ('state=0', 1),
            ('descr=description', 2)  # Should match both tasks
        ]

        for query_params, expected_count in filter_tests:
            response = self.app.get(
                f'/api/tasks?{query_params}',
                headers=self.headers
            )
            data = json.loads(response.data)

            self.assertEqual(response.status_code, 200)
            self.assertEqual(len(data['tasks']), expected_count)

        os.environ.pop('HEVELIUS_DB_NAME')

    @use_repository
    def test_tasks_performed_date_range(self, config):
        """Test tasks endpoint filtering by performed date range"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # Create test tasks with different performed dates
        # Note: We'll need to update the performed dates directly in the database
        test_task = {
            "user_id": 1,
            "scope_id": 1,
            "object": "M31",
            "ra": 0.712,
            "decl": 41.27,
            "exposure": 300.0,
            "state": 1
        }

        # Add test task
        response = self.app.post('/api/task-add',
                                 data=json.dumps(test_task),
                                 headers=self.headers)
        task_id = json.loads(response.data)['task_id']

        # Update performed date in database
        cnx = db.connect()
        db.run_query(
            cnx,
            "UPDATE tasks SET performed = %s WHERE task_id = %s",
            (datetime(2024, 1, 1, 12, 0, 0), task_id)
        )
        cnx.close()

        # Test date range filtering
        date_range_tests = [
            # Should find the task
            ('2024-01-01T00:00:00', '2024-01-02T00:00:00', 1),
            # Should not find the task
            ('2024-01-02T00:00:00', '2024-01-03T00:00:00', 0),
            # Should find the task
            ('2023-12-31T00:00:00', '2024-01-02T00:00:00', 1)
        ]

        for after, before, expected_count in date_range_tests:
            response = self.app.get(
                f'/api/tasks?performed_after={after}&performed_before={before}',
                headers=self.headers
            )
            data = json.loads(response.data)

            self.assertEqual(response.status_code, 200)
            self.assertEqual(len(data['tasks']), expected_count)

        os.environ.pop('HEVELIUS_DB_NAME')

    def test_tasks_invalid_sort_field(self):
        """Test tasks endpoint with invalid sort field"""
        response = self.app.get(
            '/api/tasks?sort_by=invalid_field',
            headers=self.headers
        )
        self.assertEqual(response.status_code, 422)  # Unprocessable Entity

    def test_tasks_invalid_sort_order(self):
        """Test tasks endpoint with invalid sort order"""
        response = self.app.get(
            '/api/tasks?sort_order=invalid',
            headers=self.headers
        )
        self.assertEqual(response.status_code, 422)  # Unprocessable Entity

    def test_tasks_no_auth(self):
        """Test tasks endpoint without authentication"""
        response = self.app.get('/api/tasks')
        self.assertEqual(response.status_code, 401)  # Unauthorized


class TestScopes(unittest.TestCase):
    def setUp(self):
        """Set up test client before each test"""
        self.app = app.test_client()
        self.app.testing = True

        # Create a test JWT token
        with app.app_context():
            self.test_token = create_access_token(
                identity=1,  # user_id=1
                additional_claims={
                    'permissions': 1,
                    'username': 'test_user'
                }
            )
            self.headers = {
                'Authorization': f'Bearer {self.test_token}',
                'Content-Type': 'application/json'
            }

    @use_repository(load_test_data=None)
    def test_scopes_success(self, config):
        """Test successful retrieval of telescopes and their sensors"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        # First create test sensors
        cnx = db.connect()
        db.run_query(cnx, """
            INSERT INTO sensors (sensor_id, name, resx, resy, pixel_x, pixel_y, bits, width, height)
            VALUES
            (1, 'Test Sensor 1', 1024, 1024, 9.0, 9.0, 16, 9.216, 9.216),
            (2, 'Test Sensor 2', 2048, 2048, 4.5, 4.5, 16, 9.216, 9.216)
            RETURNING sensor_id""")
        cnx.close()

        # Create test telescopes
        cnx = db.connect()
        db.run_query(cnx, """
            INSERT INTO telescopes (scope_id, name, descr, min_dec, max_dec, focal, aperture,
                                  lon, lat, alt, sensor_id, active)
            VALUES
            (1, 'Test Scope 1', 'Test Description 1', -90.0, 90.0, 1000.0, 200.0,
             0.0, 0.0, 0.0, 1, true),
            (2, 'Test Scope 2', 'Test Description 2', -45.0, 45.0, 2000.0, 400.0,
             0.0, 0.0, 0.0, 2, true),
            (3, 'Test Scope 3', 'Test Description 3', -30.0, 30.0, 3000.0, 600.0,
             0.0, 0.0, 0.0, NULL, false)
            RETURNING scope_id""")
        cnx.close()

        # Get scopes
        response = self.app.get('/api/scopes', headers=self.headers)
        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertIn('telescopes', data)
        telescopes = data['telescopes']

        # Should find all 3 telescopes
        self.assertEqual(len(telescopes), 3)

        # Verify first telescope with sensor
        scope1 = next(t for t in telescopes if t['scope_id'] == 1)
        self.assertEqual(scope1['name'], 'Test Scope 1')
        self.assertEqual(scope1['descr'], 'Test Description 1')
        self.assertEqual(scope1['min_dec'], -90.0)
        self.assertEqual(scope1['max_dec'], 90.0)
        self.assertEqual(scope1['focal'], 1000.0)
        self.assertEqual(scope1['aperture'], 200.0)
        self.assertEqual(scope1['lon'], 0.0)
        self.assertEqual(scope1['lat'], 0.0)
        self.assertEqual(scope1['alt'], 0.0)
        self.assertTrue(scope1['active'])
        self.assertIsNotNone(scope1['sensor'])
        self.assertEqual(scope1['sensor']['sensor_id'], 1)
        self.assertEqual(scope1['sensor']['name'], 'Test Sensor 1')
        self.assertEqual(scope1['sensor']['resx'], 1024)
        self.assertEqual(scope1['sensor']['resy'], 1024)
        self.assertEqual(scope1['sensor']['pixel_x'], 9.0)
        self.assertEqual(scope1['sensor']['pixel_y'], 9.0)
        self.assertEqual(scope1['sensor']['bits'], 16)
        self.assertEqual(scope1['sensor']['width'], 9.216)
        self.assertEqual(scope1['sensor']['height'], 9.216)

        # Verify second telescope with different sensor
        scope2 = next(t for t in telescopes if t['scope_id'] == 2)
        self.assertEqual(scope2['name'], 'Test Scope 2')
        self.assertEqual(scope2['sensor']['sensor_id'], 2)
        self.assertEqual(scope2['sensor']['name'], 'Test Sensor 2')
        self.assertEqual(scope2['sensor']['resx'], 2048)
        self.assertEqual(scope2['sensor']['resy'], 2048)
        self.assertEqual(scope2['sensor']['pixel_x'], 4.5)
        self.assertEqual(scope2['sensor']['pixel_y'], 4.5)

        # Verify third telescope without sensor
        scope3 = next(t for t in telescopes if t['scope_id'] == 3)
        self.assertEqual(scope3['name'], 'Test Scope 3')
        self.assertIsNone(scope3['sensor'])
        self.assertFalse(scope3['active'])

        os.environ.pop('HEVELIUS_DB_NAME')

    def test_scopes_no_auth(self):
        """Test scopes endpoint without authentication"""
        response = self.app.get('/api/scopes')
        self.assertEqual(response.status_code, 401)  # Unauthorized

    @use_repository(load_test_data=None)
    def test_scopes_empty(self, config):
        """Test scopes endpoint with empty database"""
        os.environ['HEVELIUS_DB_NAME'] = config['database']

        response = self.app.get('/api/scopes', headers=self.headers)
        data = json.loads(response.data)

        self.assertEqual(response.status_code, 200)
        self.assertIn('telescopes', data)
        self.assertEqual(len(data['telescopes']), 0)

        os.environ.pop('HEVELIUS_DB_NAME')


if __name__ == '__main__':
    unittest.main()
