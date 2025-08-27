"""API blueprint initialization."""

from flask import Blueprint

bp = Blueprint('api', __name__)

# Only import routes - remove analysis and export imports
from app.api import routes