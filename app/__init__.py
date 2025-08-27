"""
KinoPlex - Next-Generation Phosphoproteomics Platform
Main application factory and configuration
"""

from flask import Flask, g
from flask_cors import CORS
import logging
import os
import sqlite3


def close_db(e=None):
    """Close database connection."""
    db = g.pop('db', None)
    if db is not None:
        db.close()


def create_app(config_name=None):
    """Application factory pattern."""
    app = Flask(__name__)

    # Basic configuration
    app.config['SECRET_KEY'] = 'kinoplex-dev-secret-key-2025'
    app.config['DATABASE_PATH'] = os.environ.get('DATABASE_PATH', 'kinoplex.db')

    # Initialize extensions
    CORS(app)

    # Register teardown handler
    app.teardown_appcontext(close_db)

    # Register blueprints
    from app.main import bp as main_bp
    app.register_blueprint(main_bp)

    from app.api import bp as api_bp
    app.register_blueprint(api_bp, url_prefix='/api')

    # Configure basic logging
    if not app.debug:
        logging.basicConfig(level=logging.INFO)
        app.logger.setLevel(logging.INFO)
        app.logger.info('KinoPlex startup')

    return app