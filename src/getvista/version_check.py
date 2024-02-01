import importlib.metadata
from packaging import version

# Check for package updates
def check_for_updates():
    try:
        available_version = importlib.metadata.version('getvista')
        current_version = '0.6.2'  # Your current version from setup.py

        if version.parse(available_version) > version.parse(current_version):
            print(f"Update available: {available_version}. You are using {current_version}. Run 'pip install getvista --upgrade' to get the latest version.")
    except importlib.metadata.PackageNotFoundError:
        pass  # Package is not installed, so version check is not possible
