# Dique en Talud — Diseño del manto principal

Aplicación web para el diseño de diques en talud con manto principal de bloques de hormigón o escollera natural.

**Autores:** R. Bergillos, M. Ortega, M. Clavero — GDFA · IISTA · UGR · 2018-2025  
**Versión web:** M. Clavero Gilabert — 2026

---

## Despliegue en Streamlit Community Cloud (gratuito)

### Paso 1 — Crear el repositorio en GitHub

1. Ve a [github.com](https://github.com) → **New repository**.
2. Dale un nombre (p. ej. `dique-talud-app`), márcalo como **Public** y crea el repositorio.
3. Sube los tres archivos:
   - `app.py`
   - `requirements.txt`
   - `README.md`

   Puedes hacerlo desde la web de GitHub (botón *Add file → Upload files*) o con Git:

   ```bash
   git init
   git add app.py requirements.txt README.md
   git commit -m "Initial commit — Streamlit app dique en talud"
   git branch -M main
   git remote add origin https://github.com/TU_USUARIO/dique-talud-app.git
   git push -u origin main
   ```

### Paso 2 — Desplegar en Streamlit Community Cloud

1. Ve a [share.streamlit.io](https://share.streamlit.io) e inicia sesión con tu cuenta de GitHub.
2. Haz clic en **New app**.
3. Selecciona:
   - **Repository:** `TU_USUARIO/dique-talud-app`
   - **Branch:** `main`
   - **Main file path:** `app.py`
4. Pulsa **Deploy!**

La app tarda 1–2 minutos en arrancar la primera vez. Cuando esté lista obtendrás una URL pública del tipo:

```
https://tu-usuario-dique-talud-app-app-xxxx.streamlit.app
```

Esa URL la puedes compartir con quien quieras sin que necesiten instalar nada.

### Actualizar la app

Cada vez que hagas `git push` al repositorio, Streamlit redeployará la app automáticamente.

---

## Ejecución local

```bash
pip install -r requirements.txt
streamlit run app.py
```

La app abrirá en `http://localhost:8501`.

---

## Estructura del repositorio

```
dique-talud-app/
├── app.py            # App Streamlit (cálculos + interfaz)
├── requirements.txt  # Dependencias Python
└── README.md         # Este archivo
```
