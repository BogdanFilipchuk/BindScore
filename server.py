# server.py
from flask import Flask, request, jsonify
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt
import jwt
import datetime
import os

app = Flask(__name__)
CORS(app)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///users.db"
app.config["SECRET_KEY"] = "your-secret-key-change-this"

db      = SQLAlchemy(app)
bcrypt  = Bcrypt(app)


# ── Models ───────────────────────────────────────────────
class User(db.Model):
    id       = db.Column(db.Integer, primary_key=True)
    email    = db.Column(db.String(120), unique=True, nullable=False)
    password = db.Column(db.String(200), nullable=False)
    data     = db.relationship("UserData", backref="user", lazy=True)

class UserData(db.Model):
    id         = db.Column(db.Integer, primary_key=True)
    user_id    = db.Column(db.Integer, db.ForeignKey("user.id"), nullable=False)
    content    = db.Column(db.Text, nullable=False)
    source     = db.Column(db.String(50))   # "extension" or "desktop"
    created_at = db.Column(db.DateTime, default=datetime.datetime.utcnow)


# ── Auth helpers ─────────────────────────────────────────
def make_token(user_id):
    payload = {
        "user_id": user_id,
        "exp":     datetime.datetime.utcnow() + datetime.timedelta(days=7)
    }
    return jwt.encode(payload, app.config["SECRET_KEY"], algorithm="HS256")

def get_current_user():
    token = request.headers.get("Authorization", "").replace("Bearer ", "")
    if not token:
        return None
    try:
        payload = jwt.decode(token, app.config["SECRET_KEY"], algorithms=["HS256"])
        return User.query.get(payload["user_id"])
    except jwt.ExpiredSignatureError:
        return None
    except jwt.InvalidTokenError:
        return None


# ── Auth routes ──────────────────────────────────────────
@app.route("/register", methods=["POST"])
def register():
    body  = request.get_json()
    email = body.get("email")
    pw    = body.get("password")

    if User.query.filter_by(email=email).first():
        return jsonify({"error": "Email already registered"}), 400

    hashed = bcrypt.generate_password_hash(pw).decode("utf-8")
    user   = User(email=email, password=hashed)
    db.session.add(user)
    db.session.commit()
    return jsonify({"token": make_token(user.id)})

@app.route("/login", methods=["POST"])
def login():
    body  = request.get_json()
    email = body.get("email")
    pw    = body.get("password")
    user  = User.query.filter_by(email=email).first()

    if not user or not bcrypt.check_password_hash(user.password, pw):
        return jsonify({"error": "Invalid credentials"}), 401

    return jsonify({"token": make_token(user.id)})


# ── Data routes ──────────────────────────────────────────
@app.route("/data", methods=["POST"])
def save_data():
    """Extension or desktop app sends data here."""
    user = get_current_user()
    if not user:
        return jsonify({"error": "Unauthorized"}), 401

    body   = request.get_json()
    entry  = UserData(
        user_id = user.id,
        content = str(body.get("content")),
        source  = body.get("source", "unknown")
    )
    db.session.add(entry)
    db.session.commit()
    return jsonify({"status": "saved", "id": entry.id})

@app.route("/data", methods=["GET"])
def get_data():
    """Extension or desktop app fetches user's data."""
    user = get_current_user()
    if not user:
        return jsonify({"error": "Unauthorized"}), 401

    entries = UserData.query.filter_by(user_id=user.id)\
                            .order_by(UserData.created_at.desc()).all()
    return jsonify([{
        "id":      e.id,
        "content": e.content,
        "source":  e.source,
        "date":    e.created_at.isoformat()
    } for e in entries])


if __name__ == "__main__":
    with app.app_context():
        db.create_all()
    app.run(debug=True)