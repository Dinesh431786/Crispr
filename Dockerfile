FROM python:3.11-slim

WORKDIR /app
COPY crispr_app/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY crispr_app/ /app/

EXPOSE 8000
# Set CRISPR_CORS_ORIGINS to a comma-separated allowlist to lock down CORS.
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
