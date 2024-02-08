import argparse
import configparser
import os

def update_email(new_email):
    config = configparser.ConfigParser()
    # Get the directory path of the current script
    current_dir = os.path.dirname(__file__)
    # Specify the path to config.ini relative to the script's directory
    config_file_path = os.path.join(current_dir, 'config.ini')
    config.read(config_file_path)
    email_address = config.get('User', 'email', fallback=None)
    if email_address:
        print(f"Previously saved email address: {email_address}")
    else:
        print("No email address was previously saved.")
    config.set('User', 'email', new_email)
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)
    print(f"Email address updated successfully to {new_email}.")

def main():
    parser = argparse.ArgumentParser(description='Manage email address for GenBank Entrez queries. Only stored locally and sent with queries to NCBI, nowhere else.')
    parser.add_argument('-check', action='store_true', help='Check the current saved email address')
    parser.add_argument('-update', action='store_true', help='Update the email address')
    args = parser.parse_args()

    if args.check:
        config = configparser.ConfigParser()
        # Get the directory path of the current script
        current_dir = os.path.dirname(__file__)
        # Specify the path to config.ini relative to the script's directory
        config_file_path = os.path.join(current_dir, 'config.ini')
        config.read(config_file_path)
        email_address = config.get('User', 'email', fallback=None)
        if email_address:
            print(f"Saved email address: {email_address}")
        else:
            print("No email address saved.")

    if args.update:
        new_email = input("Enter your new email address: ")
        update_email(new_email)

if __name__ == "__main__":
    main()
