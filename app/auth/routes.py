# app/auth/routes.py
"""Authentication routes (placeholder for future implementation)."""

from flask import render_template, redirect, url_for, flash, request
from app.auth import bp


@bp.route('/login')
def login():
    """User login (placeholder)."""
    return render_template('auth/login.html', title='Login')


@bp.route('/register')
def register():
    """User registration (placeholder)."""
    return render_template('auth/register.html', title='Register')


@bp.route('/logout')
def logout():
    """User logout (placeholder)."""
    flash('You have been logged out.')
    return redirect(url_for('main.index'))