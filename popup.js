const API = "https://your-server.com";  // or http://localhost:5000 for testing

// ── Save token locally in extension ──────────────────────
async function getToken() {
  return new Promise(resolve => {
    chrome.storage.local.get("token", r => resolve(r.token));
  });
}
async function saveToken(token) {
  return new Promise(resolve => {
    chrome.storage.local.set({ token }, resolve);
  });
}

// ── Login ─────────────────────────────────────────────────
async function login(email, password) {
  const res = await fetch(`${API}/login`, {
    method:  "POST",
    headers: { "Content-Type": "application/json" },
    body:    JSON.stringify({ email, password })
  });
  const data = await res.json();
  if (data.token) {
    await saveToken(data.token);
    showMain();
  } else {
    document.getElementById("error").textContent = data.error;
  }
}

// ── Send current page to server ───────────────────────────
async function sendPage() {
  const token = await getToken();
  const [tab] = await chrome.tabs.query({ active: true, currentWindow: true });

  const res = await fetch(`${API}/data`, {
    method:  "POST",
    headers: {
      "Content-Type":  "application/json",
      "Authorization": `Bearer ${token}`
    },
    body: JSON.stringify({
      content: { url: tab.url, title: tab.title },
      source:  "extension"
    })
  });

  const result = await res.json();
  document.getElementById("status").textContent =
    result.status === "saved" ? "✓ Saved!" : "Error: " + result.error;
}

// ── On load: check if already logged in ───────────────────
window.onload = async () => {
  const token = await getToken();
  token ? showMain() : showLogin();
};

function showLogin() {
  document.getElementById("login-view").style.display = "block";
  document.getElementById("main-view").style.display  = "none";
}
function showMain() {
  document.getElementById("login-view").style.display = "none";
  document.getElementById("main-view").style.display  = "block";
}

document.getElementById("btn-login").addEventListener("click", () => {
  login(
    document.getElementById("email").value,
    document.getElementById("password").value
  );
});
document.getElementById("btn-send").addEventListener("click", sendPage);