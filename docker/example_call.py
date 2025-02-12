import json
import threading

import requests
from docker import grb_webservice
from docker.grb_webservice import start_server

# #TEST GET
# print("GET request")
# url = base_url
# x = requests.get(url)
# print(x.text)
#
# #TEST POST
# print("POST request")
# url = base_url +'/actualiser'
# with open('body_example.json', "r") as f:
#     request_body = json.load(f)
# x = requests.post(url, json = request_body)
# print(x.text)

# print("start webservice")
# grb_webservice.start_server()
# print("started webservice")


base_url = 'http://localhost:7999'
# Start the web service in a separate thread
thread = threading.Thread(target=start_server)
thread.start()

# Wait a moment for the server to start
import time
time.sleep(2)

# Make a call to the web service
url = base_url +'/actualiser'
with open('body_example.json', "r") as f:
    request_body = json.load(f)
response = requests.post(url, json=request_body)
print(response.json())