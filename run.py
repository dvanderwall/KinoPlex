#!/usr/bin/env python3
"""
KinoPlex Application Runner - Working Version
Place this at your project root and run: python run.py
"""

import os
import sys
import sqlite3
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def find_and_test_database():
    """Find and test the database connection."""
    possible_paths = [
        project_root / 'kinoplex.db',
        Path('/Users/davidvanderwall/Desktop/KinoPlex_App/kinoplex.db'),
        project_root / 'app' / 'kinoplex.db',
        Path('./kinoplex.db')
    ]

    print("Searching for kinoplex.db...")
    for path in possible_paths:
        print(f"  Checking: {path}")
        if path.exists():
            print(f"  ✓ Found: {path}")

            # Test the database
            try:
                conn = sqlite3.connect(str(path))
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM phosphosites LIMIT 1")
                count = cursor.fetchone()[0]
                conn.close()

                print(f"  ✓ Database test successful: {count:,} phosphosites")
                return str(path)

            except Exception as e:
                print(f"  ✗ Database test failed: {e}")
                continue
        else:
            print(f"  ✗ Not found")

    return None


def main():
    """Main application entry point."""
    print("🧬 KinoPlex - Next-Generation Phosphoproteomics Platform")
    print("=" * 60)

    # Find database
    db_path = find_and_test_database()
    if not db_path:
        print("\n❌ ERROR: Could not find or connect to kinoplex.db")
        print("\nPlease ensure kinoplex.db is in one of these locations:")
        print("  - Project root directory")
        print("  - /Users/davidvanderwall/Desktop/KinoPlex_App/")
        print("  - app/ subdirectory")
        sys.exit(1)

    # Set environment variable
    os.environ['DATABASE_PATH'] = db_path
    print(f"\n✅ Database configured: {db_path}")

    # Import and create app
    try:
        print("🔧 Initializing Flask application...")
        from app import create_app

        app = create_app()
        print("✅ Application created successfully!")

        # Test a simple route
        with app.test_client() as client:
            response = client.get('/api/stats/overview')
            if response.status_code == 200:
                print("✅ API endpoints working!")
            else:
                print(f"⚠️  API test returned status: {response.status_code}")

    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("\nPlease check:")
        print("  - All required files are in place")
        print("  - You're running from the project root directory")
        print("  - Flask is installed: pip install flask flask-cors")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Application error: {e}")
        sys.exit(1)

    print(f"""
🚀 Starting KinoPlex Web Server
📱 Open your browser to: http://localhost:5000
🛑 Press Ctrl+C to stop the server
{"=" * 60}
""")

    try:
        app.run(
            host='0.0.0.0',
            port=5000,
            debug=True,
            use_reloader=False  # Disable reloader to avoid import issues
        )
    except KeyboardInterrupt:
        print("\n🛑 Server stopped by user")
    except Exception as e:
        print(f"\n❌ Server error: {e}")


if __name__ == '__main__':
    main()