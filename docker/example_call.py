import json
import threading

import requests

from grb_webservice import start_server

base_url = "http://localhost:7999"

# Start the web service in a separate thread
thread = threading.Thread(target=start_server)
thread.start()

# Wait a moment for the server to start
import time

time.sleep(2)

# Make a call to the web service
url = base_url + "/actualiser"

with open("body.json", "r") as f:
    request_body = json.load(f)
response = requests.post(url, json=request_body)
print(response.json())
