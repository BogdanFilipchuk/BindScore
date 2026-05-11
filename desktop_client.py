# desktop_client.py
import requests

API = "https://your-server.com"  # same server

def login(email, password):
    res = requests.post(f"{API}/login",
                        json={"email": email, "password": password})
    return res.json().get("token")

def save_data(token, content):
    res = requests.post(f"{API}/data",
                        headers={"Authorization": f"Bearer {token}"},
                        json={"content": content, "source": "desktop"})
    return res.json()

def get_data(token):
    res = requests.get(f"{API}/data",
                       headers={"Authorization": f"Bearer {token}"})
    return res.json()

# Example usage
token = login("user@email.com", "password123")
save_data(token, {"pdb": "1BRS", "dG": -9.5})
entries = get_data(token)
for e in entries:
    print(e["source"], "→", e["content"])