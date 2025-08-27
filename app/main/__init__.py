# app/main/__init__.py
"""Main blueprint for web interface routes."""

from flask import Blueprint

bp = Blueprint('main', __name__)

from app.main import routes