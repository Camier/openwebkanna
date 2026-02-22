import os

c = get_config()


def env_bool(name: str, default: bool) -> bool:
    value = os.getenv(name, str(default)).strip().lower()
    return value in {"1", "true", "yes", "on"}


def env_int(name: str, default: int) -> int:
    try:
        return int(os.getenv(name, str(default)))
    except ValueError:
        return default


# Keep server headless and avoid extra local attack surface.
c.ServerApp.open_browser = False
c.ServerApp.ip = "0.0.0.0"
c.ServerApp.port = env_int("JUPYTER_INTERNAL_PORT", 8889)
c.ServerApp.allow_remote_access = env_bool("JUPYTER_ALLOW_REMOTE_ACCESS", False)
c.ServerApp.local_hostnames = [
    item.strip()
    for item in os.getenv(
        "JUPYTER_LOCAL_HOSTNAMES",
        "localhost,jupyter,openwebui_jupyter",
    ).split(",")
    if item.strip()
]
c.ServerApp.terminals_enabled = env_bool("JUPYTER_TERMINALS_ENABLED", False)
c.ServerApp.disable_check_xsrf = False

jupyter_token = os.getenv("JUPYTER_TOKEN", "").strip()
if not jupyter_token:
    raise RuntimeError("JUPYTER_TOKEN is required for a secure Jupyter server.")
c.ServerApp.token = jupyter_token
c.ServerApp.allow_password_change = env_bool("JUPYTER_ALLOW_PASSWORD_CHANGE", False)

web_ui_port = os.getenv("WEBUI_PORT", "3000")
default_origin_pat = rf"https?://(localhost|127\.0\.0\.1)(:{web_ui_port})?"

if env_bool("JUPYTER_ALLOW_CORS", True):
    c.ServerApp.allow_origin = os.getenv("JUPYTER_ALLOW_ORIGIN", "")
    c.ServerApp.allow_origin_pat = os.getenv(
        "JUPYTER_ALLOW_ORIGIN_PAT",
        default_origin_pat,
    )
    c.ServerApp.allow_credentials = env_bool("JUPYTER_ALLOW_CREDENTIALS", False)
else:
    c.ServerApp.allow_origin = ""
    c.ServerApp.allow_origin_pat = ""
    c.ServerApp.allow_credentials = False

# Cull idle kernels to bound resource usage.
c.MappingKernelManager.cull_idle_timeout = env_int("JUPYTER_CULL_IDLE_TIMEOUT", 900)
c.MappingKernelManager.cull_interval = env_int("JUPYTER_CULL_INTERVAL", 120)
c.MappingKernelManager.cull_connected = env_bool("JUPYTER_CULL_CONNECTED", True)
c.MappingKernelManager.cull_busy = env_bool("JUPYTER_CULL_BUSY", False)

# Optional full server shutdown when inactive (0 disables).
c.ServerApp.shutdown_no_activity_timeout = env_int("JUPYTER_SHUTDOWN_NO_ACTIVITY_TIMEOUT", 0)

c.Application.log_level = os.getenv("JUPYTER_LOG_LEVEL", "INFO")

# Prevent ad-hoc extension installs in production unless explicitly overridden.
c.LabApp.extension_manager = os.getenv("JUPYTER_LAB_EXTENSION_MANAGER", "readonly")
