"""Compatibility entrypoint.

Streamlit has been removed. Use `uvicorn main:app --reload` from this directory.
"""

from main import app

__all__ = ["app"]
