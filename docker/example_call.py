import json
import requests

url = 'http://localhost:8000/actualiser'
with open('body_example.json', "r") as f:
    request_body = json.load(f)

x = requests.post(url, json = request_body)

print(x.text)