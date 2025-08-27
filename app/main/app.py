#!/usr/bin/env python3
"""
KinoPlex Web Application - Production Runner
Main application entry point for the KinoPlex phosphoproteomics platform.
"""

from app import create_app, db
from app.database import DatabaseManager
import os
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
)
logger = logging.getLogger(__name__)


def create_production_app():
    """Create application with production configuration."""
    app = create_app()

    # Test database connection
    with app.app_context():
        try:
            # Test database connection
            db_path = app.config.get('DATABASE_PATH', 'kinoplex.db')
            if not os.path.exists(db_path):
                logger.error(f"Database file not found: {db_path}")
                return None

            # Test basic query
            test_query = "SELECT COUNT(*) as count FROM phosphosites LIMIT 1"
            result = DatabaseManager.execute_query(test_query, fetch_one=True)

            if result:
                count = result['count']
                logger.info(f"Database connected successfully. {count:,} phosphosites available.")
            else:
                logger.error("Database query failed")
                return None

        except Exception as e:
            logger.error(f"Database connection failed: {e}")
            return None

    return app


if __name__ == '__main__':
    # Set database path
    os.environ['DATABASE_PATH'] = os.path.join(os.getcwd(), 'kinoplex.db')

    # Create application
    app = create_production_app()

    if app:
        logger.info("Starting KinoPlex web application...")
        app.run(
            debug=True,
            host='0.0.0.0',
            port=5000,
            threaded=True
        )
    else:
        logger.error("Failed to initialize application. Exiting.")