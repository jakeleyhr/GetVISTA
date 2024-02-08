import importlib.metadata
import requests
from packaging import version

# Check for package updates
def check_for_updates():
    latest_version = None  # Initialize latest_version
    
    try:
        response = requests.get(f'https://pypi.org/pypi/getvista/json')
        data = response.json()
        latest_version = data['info']['version']
    except requests.RequestException:
        print('could not identify latest package version')
        pass  # Failed to fetch latest version
    
    try:
        current_version = importlib.metadata.version('getvista')
        if latest_version and version.parse(latest_version) > version.parse(current_version):
            print(f"Update available: {latest_version}. You are using {current_version}. Run 'pip install getvista --upgrade' to get the latest version.")
    except importlib.metadata.PackageNotFoundError:
        pass  # Package is not installed, so version check is not possible

if __name__ == '__main__':
    check_for_updates()