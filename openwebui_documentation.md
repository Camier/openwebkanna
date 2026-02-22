# Open WebUI Complete Documentation

Table of Contents
=================

## Getting Started
- [Quick Start](#quick-start)
  - [Important Note on User Roles and Privacy](#important-note-on-user-roles-and-privacy)
  - [Docker Installation](#docker-installation)
  - [Docker Compose Setup](#docker-compose-setup)
  - [Python Installation](#python-installation)
  - [Kubernetes Installation](#kubernetes-installation)
  - [Updating](#updating)
- [Environment Variable Configuration](#environment-variable-configuration)
  - [Important Note on PersistentConfig Environment Variables](#important-note-on-persistentconfig-environment-variables)
  - [Key Environment Variables](#key-environment-variables)
- [Advanced Topics](#advanced-topics)
  - [Understanding Open WebUI Logging](#understanding-open-webui-logging)
  - [API Keys & Monitoring](#api-keys-monitoring)
  - [Enabling HTTPS Encryption](#enabling-https-encryption)

## Features
- [Audio Features](#audio-features)
  - [Speech-to-Text](#speech-to-text)
  - [Text-to-Speech](#text-to-speech)
- [Authentication](#authentication)
  - [LDAP](#ldap)
  - [SCIM](#scim)
  - [SSO](#sso)
- [Chat Features](#chat-features)
- [Other Features](#other-features)

## Configuration
- [Configuration](#configuration)

## FAQ
- [Customization and Branding](#customization-and-branding)
- [Data Privacy and Security](#data-privacy-and-security)
- [Docker and Deployment](#docker-and-deployment)
- [HTTPS and SSL](#https-and-ssl)
- [Authentication and Sessions](#authentication-and-sessions)
- [Models and Features](#models-and-features)
- [MCP and Protocol Support](#mcp-and-protocol-support)
- [Scalability and Enterprise](#scalability-and-enterprise)
- [Release Schedule](#release-schedule)
- [License Compliance](#license-compliance)
- [Workspace](#workspace)
- [Retrieval Augmented Generation (RAG)](#retrieval-augmented-generation-rag)
- [Updating Open WebUI](#updating-open-webui)

## API Documentation
- [API Endpoints](#api-endpoints)
- [API Complete Reference](#api-complete-reference)
  - [Notable API Endpoints](#notable-api-endpoints)
  - [Ollama API Proxy Support](#ollama-api-proxy-support)
  - [Retrieval Augmented Generation (RAG) API](#retrieval-augmented-generation-rag-api)
  - [Complete API Endpoint Tree](#complete-api-endpoint-tree)
- [Swagger Documentation Links](#swagger-documentation-links)

## Troubleshooting
- [Complete Troubleshooting Guide](#complete-troubleshooting-guide)
  - [General Troubleshooting Tips](#general-troubleshooting-tips)
  - [Web Search Troubleshooting](#web-search-troubleshooting)
  - [Audio Troubleshooting](#audio-troubleshooting)
  - [RAG Troubleshooting](#rag-troubleshooting)
  - [OAuth/SSO Troubleshooting](#oauthsso-troubleshooting)
  - [Multi-Replica / High Availability / Concurrency](#multi-replica-high-availability-concurrency)
  - [Server Connectivity Issues](#server-connectivity-issues)
  - [Image Generation Troubleshooting](#image-generation-troubleshooting)
  - [Browser Compatibility](#browser-compatibility)
  - [Reset Admin Password](#reset-admin-password)
  - [Optimization, Performance & RAM Usage](#optimization-performance-ram-usage)

## Tutorials & Resources
- [Tips & Tricks](#tips-tricks)
  - [RAG Tutorial: Configuring RAG with Open WebUI Documentation](#rag-tutorial-configuring-rag-with-open-webui-documentation)
  - [Dual OAuth Configuration (Microsoft & Google)](#dual-oauth-configuration-microsoft-google)
  - [SQLite Database Overview](#sqlite-database-overview)
  - [Optimization, Performance & RAM Usage](#optimization-performance-ram-usage)
- [Integrations](#integrations)
  - [Azure OpenAI with EntraID](#azure-openai-with-entraid)
  - [Run DeepSeek R1 Dynamic 1.58-bit with Llama.cpp](#run-deepseek-r1-dynamic-158-bit-with-llamacpp)
  - [Backend-Controlled, UI-Compatible API Flow](#backend-controlled-ui-compatible-api-flow)
  - [Local LLM Setup with IPEX-LLM on Intel GPU](#local-llm-setup-with-ipex-llm-on-intel-gpu)
  - [Continue.dev VS Code Extension with Open WebUI](#continuedev-vs-code-extension-with-open-webui)
  - [Setting up with Custom CA Store](#setting-up-with-custom-ca-store)
  - [Browser Search Engine Integration](#browser-search-engine-integration)
  - [Monitor your LLM requests with Helicone](#monitor-your-llm-requests-with-helicone)
  - [Monitoring and Debugging with Langfuse](#monitoring-and-debugging-with-langfuse)
  - [LibreTranslate Integration](#libretranslate-integration)
  - [Redis Websocket Support](#redis-websocket-support)
  - [Integrate with Amazon Bedrock](#integrate-with-amazon-bedrock)
  - [Integrate with MiniMax M2.1](#integrate-with-minimax-m21)
  - [Integrate with OneDrive & SharePoint](#integrate-with-onedrive-sharepoint)
  - [Okta OIDC SSO Integration](#okta-oidc-sso-integration)
  - [Notion MCP Integration](#notion-mcp-integration)
  - [Jupyter Notebook Integration](#jupyter-notebook-integration)
  - [Firefox AI Chatbot Sidebar](#firefox-ai-chatbot-sidebar)
  - [Azure AD Domain Services (LDAPS) Integration](#azure-ad-domain-services-ldaps-integration)
  - [Iterm2 AI Integration](#iterm2-ai-integration)
- [Maintenance](#maintenance)
  - [Exporting and Importing Database](#exporting-and-importing-database)
  - [Switching to S3 Storage](#switching-to-s3-storage)
  - [Backups](#backups)
- [HTTPS Setup](#https-setup)
- [Offline Mode](#offline-mode)

## About & Community
- [Mission & Vision](#mission-vision)
- [Team](#team)
- [Licensing](#licensing)
- [Contributing](#contributing)
  - [How to Contribute](#how-to-contribute)
  - [Development Setup](#development-setup)
- [Enterprise](#enterprise)
- [Security Policy](#security-policy)
  - [Vulnerability Reporting](#vulnerability-reporting)
- [Community Resources](#community-resources)
  - [Discord Server](#discord-server)
  - [GitHub Discussions](#github-discussions)
  - [Social Media](#social-media)
- [Sponsorships](#sponsorships)
- [Roadmap](#roadmap)
  - [Interface Improvements](#interface-improvements)
  - [Information Retrieval](#information-retrieval)
  - [Community](#community)---

## Getting Started

### Quick Start

#### Important Note on User Roles and Privacy:

- **Admin Creation**: The first account created on Open WebUI gains **Administrator privileges**, controlling user management and system settings.
- **User Registrations**: Subsequent sign-ups start with **Pending** status, requiring Administrator approval for access.
- **Privacy and Data Security**: **All your data**, including login details, is **locally stored** on your device. Open WebUI ensures **strict confidentiality** and **no external requests** for enhanced privacy and security.
  - **All models are private by default.** Models must be explicitly shared via groups or by being made public. If a model is assigned to a group, only members of that group can see it. If a model is made public, anyone on the instance can see it.

Choose your preferred installation method:

- **Docker**: Officially supported and recommended for most users
- **Python**: Suitable for low-resource environments or those wanting a manual setup
- **Kubernetes**: Ideal for enterprise deployments that require scaling and orchestration

#### Docker Installation

**Step 1: Pull the Open WebUI Image**

```bash
docker pull ghcr.io/open-webui/open-webui:v0.8.3
```

**Slim Image Variant**

For environments with limited storage or bandwidth:

```bash
docker pull ghcr.io/open-webui/open-webui:main-slim
```

**Specific release version** (recommended for production):

```bash
docker pull ghcr.io/open-webui/open-webui:v0.8.3
```

**Step 2: Run the Container**

```bash
docker run -d -p 3000:8080 -v open-webui:/app/backend/data --name open-webui ghcr.io/open-webui/open-webui:v0.8.3
```

**Using GPU Support**

```bash
docker run -d -p 3000:8080 --gpus all -v open-webui:/app/backend/data --name open-webui ghcr.io/open-webui/open-webui:cuda
```

**Single-User Mode** (Disabling Login)

```bash
docker run -d -p 3000:8080 -e WEBUI_AUTH=False -v open-webui:/app/backend/data --name open-webui ghcr.io/open-webui/open-webui:v0.8.3
```

**Access the WebUI**

After the container is running, access Open WebUI at: `http://localhost:3000`

#### Docker Compose Setup

**Example `docker-compose.yml`**

```yaml
services:
  openwebui:
    image: ghcr.io/open-webui/open-webui:v0.8.3
    ports:
      - "3000:8080"
    volumes:
      - open-webui:/app/backend/data

volumes:
  open-webui:
```

**Starting the Services**

```bash
docker compose up -d
```

#### Python Installation

**Using `uv` runtime manager**

```bash
# macOS/Linux
DATA_DIR=~/.open-webui uvx --python 3.11 open-webui@latest serve

# Windows (PowerShell)
$env:DATA_DIR="C:\open-webui\data"; uvx --python 3.11 open-webui@latest serve
```

**Using Conda**

```bash
conda create -n open-webui python=3.11
conda activate open-webui
pip install open-webui
open-webui serve
```

**Using Virtual Environments**

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install open-webui
open-webui serve
```

#### Kubernetes Installation

**Prerequisites**

- Kubernetes cluster is set up
- Helm is installed

**Helm Steps**

```bash
helm repo add open-webui `https://open-webui.github.io/helm-charts`
helm repo update
helm install openwebui open-webui/open-webui
kubectl get pods
```

**Access the WebUI**

You can access Open WebUI by port-forwarding or configuring an Ingress.

#### Updating

**Option 1: Using Watchtower**

```bash
docker run --rm --volume /var/run/docker.sock:/var/run/docker.sock nickfedor/watchtower --run-once open-webui
```

**Option 2: Manual Update**

```bash
docker rm -f open-webui
docker pull ghcr.io/open-webui/open-webui:v0.8.3
docker run -d -p 3000:8080 -v open-webui:/app/backend/data --name open-webui ghcr.io/open-webui/open-webui:v0.8.3
```

### Environment Variable Configuration

Open WebUI provides a large range of environment variables that allow you to customize and configure various aspects of the application.

#### Important Note on `PersistentConfig` Environment Variables

When launching Open WebUI for the first time, all environment variables are treated equally. However, for environment variables marked as `PersistentConfig`, their values are persisted and stored internally.

After the initial launch, if you restart the container, `PersistentConfig` environment variables will no longer use the external environment variable values. Instead, they will use the internally stored values.

To disable this behavior and force Open WebUI to always use your environment variables, set `ENABLE_PERSISTENT_CONFIG` to `False`.

#### Key Environment Variables

**General Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `WEBUI_URL` | str | `http://localhost:3000` | Specifies the URL where your Open WebUI installation is reachable |
| `ENABLE_SIGNUP` | bool | `True` | Toggles user account creation |
| `ENABLE_LOGIN_FORM` | bool | `True` | Toggles email, password sign-in elements |
| `ENABLE_PASSWORD_AUTH` | bool | `True` | Allows both password and SSO authentication methods |
| `DEFAULT_LOCALE` | str | `en` | Sets the default locale for the application |
| `DEFAULT_MODELS` | str | Empty string | Sets a default Language Model |
| `DEFAULT_USER_ROLE` | str | `pending` | Sets the default role assigned to new users |
| `ENABLE_CHANNELS` | bool | `False` | Enables or disables channel support |
| `ENABLE_FOLDERS` | bool | `True` | Enables or disables the folders feature |
| `ENABLE_NOTES` | bool | `True` | Enables or disables the notes feature |
| `ENABLE_MEMORIES` | bool | `True` | Enables or disables the memory feature |
| `WEBUI_NAME` | str | `Open WebUI` | Sets the main WebUI name |
| `PORT` | int | `8080` | Sets the port to run Open WebUI from |

**Security Variables**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `WEBUI_AUTH` | bool | `True` | Enables or disables authentication |
| `ENABLE_PASSWORD_VALIDATION` | bool | `False` | Enables password complexity validation |
| `PASSWORD_VALIDATION_REGEX_PATTERN` | str | `^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)(?=.*[^\w\s]).{8,}$` | Regular expression for password validation |
| `WEBUI_SECRET_KEY` | str | `t0p-s3cr3t` | Override for JWT and encryption of sensitive data |
| `ENABLE_VERSION_UPDATE_CHECK` | bool | `True` | Enables automatic update checks |
| `OFFLINE_MODE` | bool | `False` | Disables network connections for update checks |
| `CORS_ALLOW_ORIGIN` | str | `*` | Sets allowed origins for CORS |

**Database Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `DATABASE_URL` | str | `sqlite:///${DATA_DIR}/webui.db` | Specifies the database connection URL |
| `DATABASE_POOL_SIZE` | int | `None` | Specifies the pooling strategy and size |
| `DATABASE_POOL_MAX_OVERFLOW` | int | `0` | Specifies the database pool max overflow |
| `ENABLE_DB_MIGRATIONS` | bool | `True` | Controls whether database migrations run automatically |

**Redis Settings** (Required for Multi-Worker Deployments)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `REDIS_URL` | str | - | Specifies the URL of the Redis instance |
| `REDIS_SENTINEL_HOSTS` | str | - | Comma-separated list of Redis Sentinels |
| `REDIS_CLUSTER` | bool | `False` | Connect to a Redis Cluster |
| `REDIS_KEY_PREFIX` | str | `open-webui` | Customizes the Redis key prefix |

**Ollama Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_OLLAMA_API` | bool | `True` | Enables the use of Ollama APIs |
| `OLLAMA_BASE_URL` | str | `http://localhost:11434` | Configures the Ollama backend URL |
| `OLLAMA_BASE_URLS` | str | - | Configures load-balanced Ollama backend hosts |

**OpenAI Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_OPENAI_API` | bool | `True` | Enables the use of OpenAI APIs |
| `OPENAI_API_BASE_URL` | str | `https://api.openai.com/v1` | Configures the OpenAI base API URL |
| `OPENAI_API_KEY` | str | - | Sets the OpenAI API key |
| `OPENAI_API_KEYS` | str | - | Supports multiple OpenAI API keys |

**Vector Database Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `VECTOR_DB` | str | `chroma` | Specifies which vector database system to use |

Options: `chroma`, `elasticsearch`, `milvus`, `opensearch`, `pgvector`, `qdrant`, `pinecone`, `s3vector`, `oracle23ai`, `weaviate`

**RAG Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `RAG_EMBEDDING_ENGINE` | str | - | Selects an embedding engine for RAG |
| `RAG_EMBEDDING_MODEL` | str | `sentence-transformers/all-MiniLM-L6-v2` | Sets a model for embeddings |
| `RAG_TOP_K` | int | `5` | Sets the default number of results for RAG |
| `CHUNK_SIZE` | int | `1500` | Sets the document chunk size for embeddings |
| `CHUNK_OVERLAP` | int | `150` | Specifies overlap between chunks |

**Web Search Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_WEB_SEARCH` | bool | `False` | Enable web search toggle |
| `WEB_SEARCH_ENGINE` | str | - | Specifies the web search engine to use |
| `WEB_SEARCH_RESULT_COUNT` | int | `3` | Maximum number of search results to crawl |

**Image Generation Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_IMAGE_GENERATION` | bool | `False` | Enables or disables image generation features |
| `IMAGE_GENERATION_ENGINE` | str | `openai` | Specifies the engine for image generation |
| `IMAGE_SIZE` | str | `512x512` | Sets the default output dimensions |

**OAuth Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_OAUTH_SIGNUP` | bool | `False` | Enables account creation when signing up via OAuth |
| `OAUTH_CLIENT_ID` | str | - | Sets the client ID for OIDC |
| `OAUTH_CLIENT_SECRET` | str | - | Sets the client secret for OIDC |
| `OPENID_PROVIDER_URL` | str | - | Path to the `.well-known/openid-configuration` endpoint |

**LDAP Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_LDAP` | bool | `False` | Enables or disables LDAP authentication |
| `LDAP_SERVER_HOST` | str | `localhost` | Sets the hostname of the LDAP server |
| `LDAP_SERVER_PORT` | int | `389` | Sets the port number of the LDAP server |
| `LDAP_APP_DN` | str | - | Sets the distinguished name for the LDAP application |
| `LDAP_APP_PASSWORD` | str | - | Sets the password for the LDAP application |
| `LDAP_SEARCH_BASE` | str | - | Sets the base to search for LDAP authentication |

**Audio Settings**

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `WHISPER_MODEL` | str | `base` | Sets the Whisper model for Speech-to-Text |
| `AUDIO_STT_ENGINE` | str | - | Specifies the Speech-to-Text engine |
| `AUDIO_TTS_ENGINE` | str | - | Specifies the Text-to-Speech engine |

### Advanced Topics

#### Understanding Open WebUI Logging

Logging is essential for debugging, monitoring, and understanding how Open WebUI is behaving.

**Browser Client Logging (Frontend)**

Open WebUI utilizes standard browser console logging. To access:

1. Open Developer Tools (F12 or Cmd+Opt+I on macOS)
2. Navigate to the "Console" tab

**Application Server/Backend Logging (Python)**

The backend uses Python's built-in `logging` module.

**Logging Levels**

| Level | Value | Description |
|-------|-------|-------------|
| `CRITICAL` | 50 | Severe errors that may lead to application termination |
| `ERROR` | 40 | Errors that indicate problems but app might still function |
| `WARNING` | 30 | Potential issues or unexpected situations |
| `INFO` | 20 | General informational messages |
| `DEBUG` | 10 | Detailed debugging information |

**Global Logging Level** (`GLOBAL_LOG_LEVEL`)

Sets the global logging level for all Open WebUI components.

```bash
# Example
GLOBAL_LOG_LEVEL=DEBUG
```

**App-Specific Logging Levels**

| Environment Variable | Component |
|---------------------|-----------|
| `AUDIO_LOG_LEVEL` | Audio processing |
| `COMFYUI_LOG_LEVEL` | ComfyUI Integration |
| `CONFIG_LOG_LEVEL` | Configuration Management |
| `DB_LOG_LEVEL` | Database Operations |
| `MODELS_LOG_LEVEL` | Model Management |
| `OLLAMA_LOG_LEVEL` | Ollama Backend Integration |
| `OPENAI_LOG_LEVEL` | OpenAI API Integration |
| `RAG_LOG_LEVEL` | RAG Pipeline |

#### API Keys & Monitoring

**Why Monitor?**

- Ensure Uptime
- Performance Insights
- Early Issue Detection
- Peace of Mind

**Levels of Monitoring**

1. **Basic Health Check**: Verifies if the Open WebUI service is running
2. **Model Connectivity Check**: Confirms Open WebUI can connect to configured models
3. **Model Response Testing**: Ensures models can actually process requests

**Basic Health Check Endpoint**

The `/health` endpoint is publicly accessible and returns `200 OK` when healthy.

```bash
curl https://your-open-webui-instance/health
```

#### Enabling HTTPS Encryption

HTTPS is highly recommended for security and is necessary for certain features like Voice Calls.

**Why HTTPS Matters**

- Privacy and Security
- Integrity
- Feature Compatibility (required for Voice Calls)
- Trust and User Confidence

**Choosing the Right HTTPS Solution**

- **Cloud Providers**: Load Balancers
- **Docker Container Environments**: Reverse Proxies (Nginx, Traefik, Caddy)
- **Cloudflare**: Simplified HTTPS
- **Ngrok**: Local Development HTTPS (not recommended for production)

---

## Features

### Audio Features

#### Speech-to-Text

Open WebUI supports multiple Speech-to-Text engines:

**Local Whisper (Default)**

Uses `faster-whisper` with quantization to `int8`.

Environment Variables:
- `WHISPER_MODEL`: Sets the Whisper model (default: `base`)
- `WHISPER_MODEL_DIR`: Directory to store Whisper model files
- `WHISPER_COMPUTE_TYPE`: Compute type for inference
- `WHISPER_LANGUAGE`: ISO 639-1 language code
- `WHISPER_MULTILINGUAL`: Toggle for multilingual model

**OpenAI STT**

Environment Variables:
- `AUDIO_STT_ENGINE`: Set to `openai`
- `AUDIO_STT_MODEL`: Model to use (default: `whisper-1`)
- `AUDIO_STT_OPENAI_API_BASE_URL`: Base URL for OpenAI API
- `AUDIO_STT_OPENAI_API_KEY`: API key for OpenAI

**Azure STT**

Environment Variables:
- `AUDIO_STT_ENGINE`: Set to `azure`
- `AUDIO_STT_AZURE_API_KEY`: Azure API key
- `AUDIO_STT_AZURE_REGION`: Azure region
- `AUDIO_STT_AZURE_LOCALES`: Locales for Azure STT

**Deepgram STT**

Environment Variables:
- `AUDIO_STT_ENGINE`: Set to `deepgram`
- `DEEPGRAM_API_KEY`: Deepgram API key

**Mistral STT**

Environment Variables:
- `AUDIO_STT_ENGINE`: Set to `mistral`
- `AUDIO_STT_MISTRAL_API_KEY`: Mistral API key
- `AUDIO_STT_MISTRAL_API_BASE_URL`: Mistral API base URL

#### Text-to-Speech

Open WebUI supports multiple Text-to-Speech engines:

**OpenAI TTS**

Environment Variables:
- `AUDIO_TTS_ENGINE`: Set to `openai`
- `AUDIO_TTS_MODEL`: Model to use (default: `tts-1`)
- `AUDIO_TTS_VOICE`: Voice to use (default: `alloy`)
- `AUDIO_TTS_OPENAI_API_BASE_URL`: Base URL for OpenAI API
- `AUDIO_TTS_OPENAI_API_KEY`: API key for OpenAI

**Azure TTS**

Environment Variables:
- `AUDIO_TTS_ENGINE`: Set to `azure`
- `AUDIO_TTS_AZURE_SPEECH_REGION`: Azure region
- `AUDIO_TTS_AZURE_SPEECH_OUTPUT_FORMAT`: Output format

**ElevenLabs TTS**

Environment Variables:
- `AUDIO_TTS_ENGINE`: Set to `elevenlabs`
- `ELEVENLABS_API_BASE_URL`: Base URL for ElevenLabs
- `AUDIO_TTS_API_KEY`: API key for ElevenLabs

### Authentication

#### LDAP

LDAP authentication allows users to log in using their corporate directory credentials.

Environment Variables:
- `ENABLE_LDAP`: Enable/disable LDAP authentication
- `LDAP_SERVER_LABEL`: Label for the LDAP server
- `LDAP_SERVER_HOST`: Hostname of the LDAP server
- `LDAP_SERVER_PORT`: Port number (default: `389`)
- `LDAP_ATTRIBUTE_FOR_MAIL`: Attribute to use as email
- `LDAP_ATTRIBUTE_FOR_USERNAME`: Attribute to use as username
- `LDAP_APP_DN`: Distinguished name for the LDAP application
- `LDAP_APP_PASSWORD`: Password for the LDAP application
- `LDAP_SEARCH_BASE`: Base to search for LDAP authentication
- `LDAP_USE_TLS`: Enable/disable TLS
- `LDAP_CA_CERT_FILE`: Path to CA certificate file
- `LDAP_VALIDATE_CERT`: Whether to validate the CA certificate

#### SCIM

SCIM 2.0 support for automated user and group provisioning.

Environment Variables:
- `SCIM_ENABLED`: Enable/disable SCIM 2.0 support
- `SCIM_TOKEN`: Bearer token for SCIM authentication

#### SSO

Single Sign-On support via OAuth/OIDC.

**Supported Providers**

- Google
- Microsoft
- GitHub
- Feishu
- OpenID Connect (OIDC)

**OAuth Environment Variables**

- `ENABLE_OAUTH_SIGNUP`: Enable account creation via OAuth
- `OAUTH_CLIENT_ID`: Client ID for OIDC
- `OAUTH_CLIENT_SECRET`: Client secret for OIDC
- `OPENID_PROVIDER_URL`: Path to OpenID configuration
- `OAUTH_SCOPES`: Scopes for OIDC authentication (default: `openid email profile`)
- `OAUTH_USERNAME_CLAIM`: Username claim (default: `name`)
- `OAUTH_EMAIL_CLAIM`: Email claim (default: `email`)
- `OAUTH_PICTURE_CLAIM`: Picture claim (default: `picture`)

**Provider-Specific Variables**

Google:
- `GOOGLE_CLIENT_ID`: Google OAuth client ID
- `GOOGLE_CLIENT_SECRET`: Google OAuth client secret
- `GOOGLE_OAUTH_SCOPE`: OAuth scope (default: `openid email profile`)

Microsoft:
- `MICROSOFT_CLIENT_ID`: Microsoft OAuth client ID
- `MICROSOFT_CLIENT_SECRET`: Microsoft OAuth client secret
- `MICROSOFT_CLIENT_TENANT_ID`: Tenant ID for Microsoft OAuth

GitHub:
- `GITHUB_CLIENT_ID`: GitHub OAuth client ID
- `GITHUB_CLIENT_SECRET`: GitHub OAuth client secret
- `GITHUB_CLIENT_SCOPE`: OAuth scope (default: `user:email`)

### Chat Features

- **Autocomplete**: Enable with `ENABLE_AUTOCOMPLETE_GENERATION`
- Chat Parameters: Configure advanced parameters for model responses
- Chat Sharing: Share conversations with others
- **Code Execution**: Execute code within chats
  - Artifacts: Interactive code artifacts
  - Mermaid: Diagram generation
  - Python: Python code execution
- Conversation Organization: Folders for organizing chats
- History Search: Search through conversation history
- Multi-Model Chats: Use multiple models in a single conversation
- Reasoning Models: Support for reasoning-capable models
- Temporal Awareness: Time-aware conversations
- URL Parameters: Configure chats via URL parameters

### Other Features

- Evaluation: Compare model outputs
- **Experimental Features**: Direct connections for tool servers
- **Image Generation**: Support for DALL-E, ComfyUI, AUTOMATIC1111, Gemini
- Interface Customization: Banners, webhooks
- **MCP Support**: Model Context Protocol integration
- **Memory**: Long-term memory for conversations
- Notes: Personal note-taking
- Pipelines: Custom processing pipelines
  - Filters: Filter pipeline stages
  - Pipes: Pipeline stages
  - Tutorials: Pipeline guides
  - Valves: Pipeline configuration
- **RAG**: Retrieval-Augmented Generation
- **RBAC**: Role-Based Access Control
- **Web Search**: Multiple search engine integrations
- **Workspace**: Shared models, knowledge, prompts, tools

---

## Configuration

Open WebUI can be configured through:

1. **Environment Variables**: Set before starting the application
2. **Admin Panel UI**: Change settings through the web interface
3. **PersistentConfig**: Some settings are stored in the database

**Configuration Priority**

For `PersistentConfig` variables:
1. Database-stored values (highest priority when enabled)
2. Environment variables
3. Default values

For regular variables:
1. Environment variables
2. Default values

---

## Troubleshooting

### Common Issues

- **Audio**: Issues with speech-to-text or text-to-speech
- **Compatibility**: Version compatibility issues
- **Connection Error**: Problems connecting to model providers
- **Image Generation**: Issues with image generation
- **Manual Database Migration**: Database migration problems
- **Multi-Replica**: Issues with multiple replicas
- **Password Reset**: Resetting user passwords
- **Performance**: Performance optimization
- **RAG**: Retrieval-Augmented Generation issues
- **SSO**: Single Sign-On problems
- **Web Search**: Web search integration issues

---

## Tutorials

### Deployment

Deployment tutorials for various platforms and configurations.

### HTTPS Setup

Step-by-step guides for setting up HTTPS:

- **Caddy**: Automatic HTTPS with Caddy
- **HAProxy**: HAProxy configuration
- **Nginx**: Nginx reverse proxy setup

### Integrations

Integration guides for:

- **Amazon Bedrock**
- **Azure AD DS LDAP**
- **Azure OpenAI**
- **Backend Controlled UI Compatible Flow**
- **Browser Search Engine**
- **Continue Dev**
- **Custom CA**
- **DeepSeekR1 Dynamic**
- **Entra Group Name Sync**
- **Firefox Sidebar**
- **Helicone**
- **IPEX_LLM**
- **iTerm2**
- **Jupyter**
- **Langfuse**
- **Libre Translate**
- **MCP Notion**
- **MiniMax**
- **Okta OIDC SSO**
- **OneDrive SharePoint**
- **Redis**

### Maintenance

- **Backups**: Backup strategies
- **Database**: Database management
- **S3 Storage**: S3 storage configuration

### Tips

- **Contributing Tutorial**: How to contribute
- **Dual OAuth Configuration**: Multiple OAuth providers
- **One-Click Ollama Launcher**: Quick Ollama setup
- **RAG Tutorial**: Retrieval-Augmented Generation guide
- **SQLite Database**: SQLite database tips
---

---

## About

### Mission

Open WebUI's mission to democratize AI access.

### License

Open WebUI License information.

### Roadmap

Future development plans and features.

### Team

The Open WebUI development team.

### Contributing

How to contribute to Open WebUI development.

### Sponsorships

Sponsorship opportunities.

### Enterprise

Open WebUI Enterprise offerings:

- **Architecture**: Enterprise architecture overview
- **Customers**: Enterprise customer stories
- **Customization**: Enterprise customization options
- **Integration**: Enterprise integrations
- **Partners**: Enterprise partners
- **Security**: Enterprise security features
- **Support**: Enterprise support options

---
### FAQ

Frequently Asked Questions about Open WebUI.

## API Endpoints

Open WebUI provides a comprehensive API that is compatible with the OpenAI API standard.

### Authentication

To ensure secure access to the API, authentication is required. You can authenticate your API requests using the Bearer Token mechanism. Obtain your API key from Settings > Account in the Open WebUI, or alternatively, use a JWT (JSON Web Token) for authentication.

#### Sign In (Get JWT Token)

```http
POST /api/v1/auths/signin
Content-Type: application/json
```

**Request Body:**
```json
{
  "email": "user@example.com",
  "password": "your-password"
}
```

**Response (200 OK):**
```json
{
  "token": "eyJhbGciOiJIUzI1NiIs...",
  "id": "user-uuid",
  "email": "user@example.com",
  "name": "User Name",
  "role": "user",
  "profile_image_url": "/user.png"
}
```

**Error Responses:**
- `400 Bad Request`: Missing email or password
- `401 Unauthorized`: Invalid credentials
- `500 Internal Server Error`: Server error

#### Sign Up (Create Account)

```http
POST /api/v1/auths/signup
Content-Type: application/json
```

**Request Body:**
```json
{
  "name": "New User",
  "email": "newuser@example.com",
  "password": "secure-password"
}
```

**Response (201 Created):**
Same as sign in response.

**Error Responses:**
- `400 Bad Request`: Invalid input data
- `409 Conflict`: Email already exists

#### Get Current User Info

```http
GET /api/v1/users/user/info
Authorization: Bearer <token>
```

**Response (200 OK):**
```json
{
  "id": "user-uuid",
  "email": "user@example.com",
  "name": "User Name",
  "role": "user",
  "profile_image_url": "/user.png"
}
```

#### API Key Authentication

API keys can be used instead of JWT tokens. Generate API keys in Settings > Account.

```http
GET /api/models
Authorization: Bearer sk-xxxxxxxxxxxxxxxxxxxxxxxx
```

**Requirements:**
1. Admin must enable "Enable API Keys" in Admin Panel > Settings
2. User must have `features.api_keys` permission

#### Token Expiration

JWT tokens expire based on `JWT_EXPIRES_IN` environment variable:
- `"-1"`: Never expire (not recommended for production)
- `"4w"`: 4 weeks
- `"24h"`: 24 hours
- `"30m"`: 30 minutes

When a token expires, you'll receive:
```json
{
  "detail": "Invalid token or expired token."
}
```
With HTTP 401 status code.

### Base URL

All API endpoints are prefixed with `/api` (e.g., `http://localhost:3000/api`).

### Notable API Endpoints

#### Retrieve All Models

- **Endpoint**: `GET /api/models`
- **Description**: Fetches all models created or added via Open WebUI.

#### Chat Completions

- **Endpoint**: `POST /api/chat/completions`
- **Description**: Serves as an OpenAI API compatible chat completion endpoint
- **Request Body**:
  ```json
  {
    "model": "model-id",
    "messages": [
      {"role": "user", "content": "Hello"}
    ],
    "stream": true
  }
  ```

#### Retrieval Augmented Generation (RAG)

The RAG system allows you to upload and process files for knowledge retrieval.

##### Uploading Files

- **Endpoint**: `POST /api/v1/files/`
- **Description**: Uploads a file to be processed by the RAG system
- **Request**: Multipart form data with the file

##### Checking File Processing Status

- **Endpoint**: `GET /api/v1/files/{id}/process/status`
- **Description**: Checks if a file has been processed and is ready for use
- **Response**:
  ```json
  {
    "status": "completed",
    "file_id": "file-id"
  }
  ```

##### Adding Files to Knowledge Collections

- **Endpoint**: `POST /api/v1/knowledge/{id}/file/add`
- **Description**: Adds a processed file to a knowledge collection
- **Request Body**:
  ```json
  {
    "file_id": "file-id"
  }
  ```

### Vector Database Configuration

Open WebUI supports multiple vector databases for RAG operations:

- **ChromaDB** (default)
- **Milvus**
- **PGVector** (PostgreSQL)
- **Qdrant**
- **Pinecone**
- **Elasticsearch**
- **OpenSearch**
- **S3Vector**
- **Oracle23ai**
- **Weaviate**

Configure via environment variables:

```bash
VECTOR_DB=chromadb
```

## Troubleshooting

### Connection Errors

When using Open WebUI behind a reverse proxy or HTTPS, you may encounter connection errors. Here's the required configuration:

#### Required Environment Variables

```bash
# Main URL Configuration
WEBUI_URL=https://your-open-webui-domain.com

# CORS Configuration (space-separated list)
CORS_ALLOW_ORIGIN="https://yourdomain.com;http://yourdomain.com;https://yourip;http://localhost:3000"

# Cookie Security Settings
WEBUI_SESSION_COOKIE_SECURE=true
WEBUI_AUTH_COOKIE_SECURE=true
WEBUI_SESSION_COOKIE_SAME_SITE=lax
WEBUI_AUTH_COOKIE_SAME_SITE=lax

# WebSocket Support
ENABLE_WEBSOCKET_SUPPORT=true
WEBSOCKET_MANAGER=redis
WEBSOCKET_REDIS_URL=redis://redis:6379/1
```

#### Nginx Configuration Example

```nginx
location / {
    proxy_pass http://localhost:3000;
    proxy_http_version 1.1;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;

    # WebSocket timeout settings
    proxy_connect_timeout 7d;
    proxy_send_timeout 7d;
    proxy_read_timeout 7d;
}
```

#### Common Issues

1. **WebSocket Connection Failed**
   - Ensure `ENABLE_WEBSOCKET_SUPPORT=true` is set
   - Configure Redis for WebSocket manager in multi-worker setups
   - Check reverse proxy WebSocket settings

2. **CORS Errors**
   - Add your domain to `CORS_ALLOW_ORIGIN`
   - Ensure protocol (http/https) matches
   - Include all required origins (IP, domain, localhost)

3. **HTTPS/SSL Issues**
   - Set cookie secure flags to `true`
   - Configure `SameSite` attribute correctly
   - Ensure reverse proxy handles SSL properly

### Image Generation

Open WebUI supports image generation through three backends: AUTOMATIC1111, ComfyUI, and OpenAI DALL-E.

#### AUTOMATIC1111 Setup

1. Ensure AUTOMATIC1111 is installed
2. Launch with API access: `./webui.sh --api --listen`
3. Configure in Open WebUI: Admin Panel > Settings > Images
4. Set Image Generation Engine to `Default (Automatic1111)`
5. Enter API URL: `http://<your_address>:7860/`
   - For Docker: Use `http://host.docker.internal:7860/`

Docker with preset environment variables:
```bash
docker run -d -p 3000:8080 --add-host=host.docker.internal:host-gateway \
  -e AUTOMATIC1111_BASE_URL=http://host.docker.internal:7860/ \
  -e ENABLE_IMAGE_GENERATION=True \
  -v open-webui:/app/backend/data --name open-webui --restart always \
  ghcr.io/open-webui/open-webui:v0.8.3
```

#### ComfyUI Setup

1. Download and extract ComfyUI from GitHub
2. Start ComfyUI with: `python main.py --listen`
3. For low VRAM: Use additional flags to reduce memory usage

Docker with preset environment variables:
```bash
docker run -d -p 3000:8080 --add-host=host.docker.internal:host-gateway \
  -e COMFYUI_BASE_URL=http://host.docker.internal:7860/ \
  -e ENABLE_IMAGE_GENERATION=True \
  -v open-webui:/app/backend/data --name open-webui --restart always \
  ghcr.io/open-webui/open-webui:v0.8.3
```

##### Setting Up FLUX.1 Models

1. **Model Checkpoints**: Download `FLUX.1-schnell` or `FLUX.1-dev` from HuggingFace
   - Place in both `models/checkpoints` and `models/unet` directories
2. **VAE Model**: Download `ae.safetensors` and place in `models/vae`
3. **CLIP Model**: Download `clip_l.safetensors` and place in `models/clip`
4. **T5XXL Model**: Download `t5xxl_fp16.safetensors` or `t5xxl_fp8_e4m3fn.safetensors` and place in `models/clip`

**Integration Steps**:
1. Admin Panel > Settings > Images > Choose `ComfyUI`
2. Enter API URL: `http://<address>:8188/`
3. Enable Dev Mode in ComfyUI (gear icon)
4. Export workflow in API format (`Save (API Format)`)
5. Upload `workflow_api.json` to Open WebUI
6. Map ComfyUI Workflow Nodes to imported workflow node IDs
7. Click Save

#### SwarmUI Integration

SwarmUI uses ComfyUI as backend:
1. Append `ComfyBackendDirect` to ComfyUI Base URL
2. Setup SwarmUI with LAN access
3. API URL format: `http://<address>:7801/ComfyBackendDirect`

#### DALL-E Integration

Configure OpenAI's DALL-E for image generation:

1. Obtain API key from OpenAI
2. Admin Panel > Settings > Images > Choose `Open AI (Dall-E)`
3. Enter OpenAI API key
4. Choose DALL-E model:
   - **DALL-E 2**: Supports `256x256`, `512x512`, `1024x1024`
   - **DALL-E 3**: Supports `1024x1024`, `1792x1024`, `1024x1792`

```bash
ENABLE_IMAGE_GENERATION=true
IMAGE_GENERATION_ENGINE=dalle
```

#### Using Image Generation

1. Use a text generation model to write an image generation prompt
2. Click the Picture icon after response completes
3. Image will be returned automatically in chat
4. You can also edit the LLM's response and use that as your prompt

### Memory & Personalization

The Memory system allows models to remember facts, preferences, and context across conversations. Currently in **Beta/Experimental** stage.

#### How It Works

The memory system stores information about you (e.g., preferences, location) and uses it in two ways:

1. **Manual Management**: Users can add, edit, or delete memories in Settings > Personalization > Memory
2. **Native Memory Tools (Agentic Mode)**: With Native Function Calling enabled, models can manage memory autonomously

##### Memory Tools

When using models with Native Function Calling (Agentic Mode):

- **`add_memory`**: Model proactively saves new facts learned about you
- **`search_memories`**: Model searches memory bank for relevant context (returns unique IDs, default 5 results)
- **`replace_memory_content`**: Model updates or corrects specific existing memory using its ID

##### Best Models for Memory

Autonomous memory management works best with frontier models:
- GPT-5
- Claude 4.5+
- Gemini 3+
- MiniMax M2.1

Small local models may struggle with appropriate memory selection.

#### Administrative Controls

##### Global Toggle
- **Admin UI**: Admin Panel > Settings > General > Features > Memories
- **Environment Variable**: `ENABLE_MEMORIES` (Default: `True`)

##### Granular Permissions
- **Admin UI**: Admin Panel > Users > Permissions > Features > Memories
- **Environment Variable**: `USER_PERMISSIONS_FEATURES_MEMORIES` (Default: `True`)

#### Privacy & Security

- Memories are stored locally in your Open WebUI database
- Specific to your user account
- Never shared across users
- Can be cleared at any time

### Model Context Protocol (MCP)

Open WebUI natively supports **MCP (Model Context Protocol)** starting in v0.6.31.

#### Prerequisites

**IMPORTANT**: You MUST set `WEBUI_SECRET_KEY` environment variable. Without it, OAuth-connected MCP tools will break every time you restart.

#### Quick Start

1. Open Admin Settings > External Tools
2. Click + (Add Server)
3. Set Type to **MCP (Streamable HTTP)**
4. Enter Server URL and Auth details
5. Save and restart if prompted

**Warning**: Make sure Type is set to **MCP (Streamable HTTP)**, not OpenAPI. Entering MCP-style configuration into OpenAPI will cause UI crashes.

#### MCP vs OpenAPI

**Choose OpenAPI if you want**:
- Enterprise readiness with SSO, API gateways, audit
- Operational resilience with standard HTTP verbs
- First-class observability and policy integration

**Choose MCP (Streamable HTTP) if you need**:
- Common tool protocol already used by your MCP servers
- Streamed protocol communication over HTTP

#### Configuration Best Practices

##### Authentication Modes

- **None**: For local MCP servers or internal networks
  - Default to "None" unless server requires token
  - Selecting "Bearer" without key sends empty Authorization header
- **Bearer**: Only if MCP server requires specific API token (must populate "Key" field)
- **OAuth 2.1**: For enterprise deployments with Identity Provider flows

##### Connection URLs

For Docker with MCP server on host machine:
- Use `http://host.docker.internal:<port>` instead of `localhost`

##### Function Name Filter List

- Default: Leave empty to expose all tools
- Workaround: If errors with empty list, try adding single comma (`,`)

#### Troubleshooting

##### "Failed to connect to MCP server"

**Solutions**:
1. Check authentication - switch to "None" if no token needed
2. Try adding comma to Function Name Filter List if empty

##### Infinite loading screen

**Cause**: Configured MCP server using OpenAPI connection type

**Solution**:
1. Admin Settings > External Tools
2. Disable/delete problematic tool
3. Refresh page (Ctrl+F5)
4. Re-add with correct Type: MCP (Streamable HTTP)

### Pipelines

Pipelines is a UI-agnostic OpenAI API plugin framework for extending Open WebUI with custom workflows.

#### When to Use Pipelines

**DO NOT USE PIPELINES IF**:
- You only need to add support for additional providers like Anthropic
- You need basic filters
- Use Open WebUI Functions instead (built-in, easier to configure)

**Use Pipelines for**:
- Computationally heavy tasks to offload from main instance
- Custom Python library integration
- Complex workflows requiring external processing

#### Examples

- Function Calling Pipeline
- Custom RAG Pipeline
- Message Monitoring with Langfuse
- Rate Limit Filter
- Real-Time Translation with LibreTranslate
- Toxic Message Filter
- And much more!

#### Quick Start with Docker

**Warning**: Pipelines has arbitrary code execution - don't fetch from untrusted sources.

```bash
# Run Pipelines container
docker run -d -p 9099:9099 --add-host=host.docker.internal:host-gateway \
  -v pipelines:/app/pipelines --name pipelines --restart always \
  ghcr.io/open-webui/pipelines:main

# Connect to Open WebUI
# Admin Panel > Settings > Connections > Add connection
# Set API URL to http://localhost:9099
# Set API key to 0p3n-w3bu!
```

**Note**: If Open WebUI is in Docker, replace `localhost` with `host.docker.internal`.

#### Custom Pipeline Installation

```bash
docker run -d -p 9099:9099 --add-host=host.docker.internal:host-gateway \
  -e PIPELINES_URLS="https://github.com/open-webui/pipelines/blob/main/examples/filters/detoxify_filter_pipeline.py" \
  -v pipelines:/app/pipelines --name pipelines --restart always \
  ghcr.io/open-webui/pipelines:main
```

#### Directory Structure

- `/pipelines` directory is core of setup
- All pipelines in this directory are automatically loaded
- Can change location with `PIPELINES_DIR` env variable
- Find examples in `https://github.com/open-webui/pipelines/blob/main/examples`

### Roadmap

#### Interface Enhancements

- **Packaged Single Binary Executable**: Simplified deployment across platforms
- **User Page**: Personal page with posts, followers, likes, comments for sharing configurations
- **AI Powered Notes**: Notion/Obsidian-style note-taking with AI integration
- **Advanced User Tracking**: Tools for tracking performance and managing costs
- **AI Workflow Tool**: Node-based tool for orchestrating AI systems
- **Local Text-to-Speech Integration**: Natural audio conversion
- **Wakeword Detection**: Hands-free voice commands
- **Better Code Execution**: Enhanced coding functionalities
- **Code Interpreter Function**: Multi-language code interpreter
- **Streamlined Fine-tune Support**: Automated dataset building from interface usage
- **Accessibility Enhancements**: Full accessibility for diverse abilities

#### Information Retrieval (RAG)

- **Customizable RAG Framework**: Drag-and-drop modular components
- **Advanced Integration**: Latest State-of-the-Art retrieval methods
- **Dedicated R&D**: Research collaboration with leading entities

#### Community

- **LLM Leaderboard**: Community-driven model evaluation and ranking
- **Sub-Communities**: Niche communities for specific datasets and prompt optimization
- **Community Whitepaper**: Detailed strategy for shaping AI platforms

#### Vision for the Future

The core vision is to enhance daily life by automating routine tasks, freeing time for what truly matters. As AI models become more capable, mundane responsibilities become automated, allowing collective time to be redistributed toward grand endeavors like space exploration, longevity research, or extra moments with loved ones.

### Performance Optimization

#### Dedicated Task Models

Offload non-chat tasks to lightweight models:

```bash
# External model for tasks
TASK_MODEL=gpt-5-nano
TASK_MODEL_EXTERNAL=gpt-5-nano

# Model providers for tasks
TASK_MODEL_PROVIDERS=[{"name": "OpenAI", "url": "https://openai.com"}]
```

#### Caching & Latency Optimization

```bash
# Base Models Cache
ENABLE_BASE_MODELS_CACHE=True
MODELS_CACHE_TTL=300

# Query Caching
ENABLE_QUERIES_CACHE=True

# RAG Context Optimization
RAG_SYSTEM_CONTEXT=True
```

#### Database Optimization

```bash
# PostgreSQL for production
DATABASE_URL=postgres://user:password@localhost:5432/webui

# Disable real-time chat saving for performance
ENABLE_REALTIME_CHAT_SAVE=False

# Session sharing settings
DATABASE_ENABLE_SESSION_SHARING=False

# Connection pooling
DATABASE_POOL_SIZE=15
DATABASE_POOL_MAX_OVERFLOW=20
```

#### High-Concurrency Optimization

```bash
# Streaming chunk size
CHAT_RESPONSE_STREAM_DELTA_CHUNK_SIZE=7

# Thread pool size
THREAD_POOL_SIZE=2000
```

#### Hardware Acceleration

```bash
# CUDA acceleration
CUDA_VISIBLE_DEVICES=0
ENABLE_CUDA=1

# ROCm (AMD)
ENABLE_ROCM=1

# CPU optimization
OMP_NUM_THREADS=4
```

### Web Search

#### Proxy Configuration

If you need to use a proxy for web search:

```bash
# HTTP proxy
HTTP_PROXY=http://proxy.example.com:8080
HTTPS_PROXY=http://proxy.example.com:8080

# No proxy for these hosts
NO_PROXY=localhost,127.0.0.1
```

#### SSL Verification

```bash
# Disable SSL verification (not recommended for production)
WEB_SEARCH_SSL_VERIFY=false
```

#### Search Engine Configuration

```bash
# Enable web search
ENABLE_WEB_SEARCH=true

# Default search engine
WEB_SEARCH_ENGINE=searxng

---
# Custom SearXNG instance
SEARXNG_QUERY_URL=https://searxng.example.com/search?q={query}
```

## Tutorials

### Deployment

The Deployment tutorials provide guidance on various deployment scenarios including Docker, Kubernetes, cloud platforms, and high-availability configurations.

### HTTPS Configuration

For production deployments, you should configure HTTPS. See the Troubleshooting section for detailed configuration examples using Nginx and other reverse proxies.

### Integrations

Open WebUI supports numerous integrations:

#### Web Search
Integration with various search engines including SearXNG, Google PSE, Brave, Kagi, Tavily, Perplexity, and more.

#### Image Generation
Support for DALL-E, ComfyUI, AUTOMATIC1111, and other image generation engines.

#### Database Backends
Integration with PostgreSQL, Redis, and various vector databases for RAG functionality.

#### Monitoring
Integration with Langfuse and other monitoring tools for debugging and performance tracking.

### Tips & Tricks

#### Contributing Tutorials
Community contributions for tutorials are welcome. See the Contributing section for details.

#### Open WebUI RAG Tutorial
Comprehensive guide for implementing Retrieval-Augmented Generation with Open WebUI.

#### Reduce RAM Usage
Tips and configurations for minimizing memory usage in resource-constrained environments.
---

### Migration

Guides for migrating from previous versions and from other platforms.

## FAQ

### Customization and Branding

**Q: How do I customize the logo and branding?**

A: You can customize the theme, logo, and branding with the Enterprise License, which unlocks exclusive enterprise features.

### Data Privacy and Security

**Q: Is my data being sent anywhere?**

A: No, your data is never sent anywhere unless you explicitly choose to share it or you connect an external model provider. Everything inside Open WebUI runs and is stored locally on your machine or server.

**Q: Can I use Open WebUI in outer space (e.g., Mars and beyond) or other extreme environments?**

A: Yes. Open WebUI is fully self-hosted and does not rely on persistent internet connectivity, making it suitable for environments where cloud-based systems are impractical.

**Q: Why am I asked to sign up? Where are my data being sent to?**

A: We require you to sign up to become the admin user for enhanced security. All information stays within your server and never leaves your device.

### Docker and Deployment

**Q: Why can't my Docker container connect to services on the host using `localhost`?**

A: Inside a Docker container, `localhost` refers to the container itself. Use the DNS name `host.docker.internal` instead.

**Q: How do I make my host's services accessible to Docker containers?**

A: Configure services to listen on all network interfaces using the IP address `0.0.0.0`, instead of `127.0.0.1`.

**Q: Why isn't my Open WebUI updating?**

A: You must pull the latest image, stop and remove the existing container, then start a new one:

```bash
docker pull ghcr.io/open-webui/open-webui:v0.8.3
docker stop open-webui
docker rm open-webui
docker run -d -p 3000:8080 -v open-webui:/app/backend/data --name open-webui --restart always ghcr.io/open-webui/open-webui:v0.8.3
```

**Q: Will I lose my data if I delete my container?**

A: Your data is safe ONLY if you have a Volume configured with `-v open-webui:/app/backend/data`. Without this, data is stored inside the container and deleting it results in permanent data loss.

**Q: Should I use the distro-packaged Docker or the official Docker package?**

A: We recommend using the official Docker package for the latest features, bug fixes, and support for features like `host.docker.internal`.

**Q: Is GPU support available in Docker?**

A: GPU support is officially available in Docker for Windows and Docker Engine on Linux. Docker Desktop for Linux and MacOS do not currently offer GPU support.

### HTTPS and SSL

**Q: Why doesn't Speech-to-Text (STT) and Text-to-Speech (TTS) work in my deployment?**

A: Modern browsers restrict STT and TTS to work only under secure HTTPS connections. Ensure your deployment uses HTTPS.

**Q: Why doesn't Open WebUI include built-in HTTPS support?**

A: The project leaves HTTPS implementation to users for maximum flexibility and customization for different environments.

**Q: Why can't Open WebUI start with an SSL error?**

A: This is likely due to the absence of SSL certificates or incorrect HuggingFace configuration. Set up a mirror:

```bash
docker run -d -p 3000:8080 -e HF_ENDPOINT=https://hf-mirror.com/ --add-host=host.docker.internal:host-gateway -v open-webui:/app/backend/data --name open-webui --restart always ghcr.io/open-webui/open-webui:v0.8.3
```

### Authentication and Sessions

**Q: I updated/restarted and now I'm being logged out, or getting "Error decrypting tokens"?**

A: Set a persistent `WEBUI_SECRET_KEY` in your environment variables to prevent session invalidation.

**Q: I updated/restarted and now my login isn't working anymore.**

A: Ensure your Docker container has a volume mounted for `/app/backend/data` to persist your data.

**Q: I tried to login and couldn't, made a new account and now I'm being told my account needs to be activated by an admin.**

A: The first account is automatically designated as the admin account. If you forget the admin password, see the Resetting the Admin Password guide.

### Models and Features

**Q: Why are my reasoning model's thinking blocks showing as raw text?**

A: Customize the thinking tags in the model's Advanced Parameters. See the Reasoning & Thinking Models guide.

**Q: RAG with Open WebUI is very bad or not working. Why?**

A: If using Ollama, increase the context length from the default 2048 tokens to 8192+ tokens to allow retrieved documents to contribute effectively.

**Q: I'm getting "The content provided is empty" when uploading files via the API.**

A: Content extraction happens asynchronously. Poll the file status endpoint until processing is complete:

```python
import requests
import time

def wait_for_processing(token, file_id):
    url = f'http://localhost:3000/api/v1/files/{file_id}/process/status'
    headers = {'Authorization': f'Bearer {token}'}

    while True:
        status = requests.get(url, headers=headers).json().get('status')
        if status == 'completed':
            return True
        elif status == 'failed':
            raise Exception("Processing failed")
        time.sleep(2)
```

**Q: I asked the model what it is and it gave the wrong answer.**

A: LLMs do not reliably know their own identity. Check the model selector in the interface or Admin Panel to verify which model you're using.

**Q: Why am I seeing multiple API requests when I only send one message?**

A: Open WebUI uses Task Models for background features like title generation, tag generation, and query generation. Configure a Task Model to use a smaller, cheaper model for these tasks.

### MCP and Protocol Support

**Q: Is MCP (Model Context Protocol) supported in Open WebUI?**

A: Yes, Open WebUI includes native support for MCP Streamable HTTP. For other MCP transports, use the official proxy adapter MCPO.

**Q: Why doesn't Open WebUI support [Specific Provider]'s latest API?**

A: Open WebUI is built around universal protocols like the OpenAI Chat Completions protocol to remain backend-agnostic and compatible with dozens of providers.

### Scalability and Enterprise

**Q: Is Open WebUI scalable for large organizations?**

A: Yes, Open WebUI is architected for massive scalability and is trusted in deployments supporting tens or hundreds of thousands of seats.

**Q: How can I deploy Open WebUI in a highly available, large-scale production environment?**

A: Use multiple containers behind a load balancer, external databases, enterprise authentication integration, and monitoring tools.

### Release Schedule

**Q: How often is Open WebUI updated?**

A: The project aims to ship major releases weekly, with bug fixes and minor updates delivered as needed.

### License Compliance

**Q: Where do I report non-compliant Open WebUI deployments?**

A: Email reports@openwebui.com with relevant details (screenshots, URLs, description of usage).

### Workspace

The Workspace provides a comprehensive environment for managing AI interactions and configurations.

#### Resources

- **Models** - Create and manage custom models tailored to specific purposes
- **Knowledge** - Manage knowledge bases for retrieval augmented generation
- **Prompts** - Create and organize reusable prompts

Each section is designed to give fine-grained control over your Open WebUI experience.

### Retrieval Augmented Generation (RAG)

**IMPORTANT FOR OLLAMA USERS**: Ollama defaults to a 2048-token context length. This severely limits RAG performance because retrieved data may not be used at all. **Increase context length to 8192+ tokens** in Admin Panel > Models > Settings > Advanced Parameters.

#### Overview

RAG enhances conversational capabilities by incorporating context from diverse sources including:
- Local and remote documents
- Web content
- Multimedia sources (e.g., YouTube videos)

Retrieved text is combined with a predefined RAG template and prefixed to the user's prompt.

#### Local and Remote RAG Integration

**Local Documents**:
1. Upload via Documents section in Workspace
2. Access using `#` symbol before query
3. Document icon appears above "Send a message" when selected

**Web Content**:
1. Start query with `#` followed by URL
2. Click formatted URL in box above chat
3. Document icon appears when successfully retrieved

**Tip**: Link to raw or reader-friendly versions for better results.

#### Context Length Configuration

For web content integration, ensure sufficient context length:
- **Ollama models**: Navigate to Admin Panel > Models > Settings > Advanced Parameters and increase to 8192+ (preferably 16000+) tokens
- **OpenAI/other models**: Use models with sufficient built-in context (e.g., GPT-4 Turbo with 128k tokens)

#### RAG Template Customization

Customize the RAG template from Admin Panel > Settings > Documents menu.

**Markdown Header Splitting**:
- Documents are first split by markdown headers (H1-H6)
- Preserves document structure
- Sections under same header kept together
- Further processed by standard character/token splitter

**Chunk Min Size Target**:
- Found in Admin Panel > Settings > Documents
- Intelligently merges small sections after markdown splitting
- Improves retrieval coherence
- Reduces total number of vectors in database

#### Chunking Configuration

Fine-tune how documents are split into chunks for embedding:

- **Chunk Size**: Maximum number of characters (or tokens) per chunk
- **Chunk Overlap**: Content shared between adjacent chunks to maintain context
- **Chunk Min Size Target**: Intelligently merges small pieces with neighbors

**Benefits of Chunk Min Size Target**:
- Improves RAG quality (eliminates tiny fragments)
- Reduces vector database size
- Speeds up retrieval & embedding
- Testing shows 90% reduction in chunk counts while improving accuracy

#### Merging Algorithm

The merging algorithm uses a single forward pass:
1. Start with first chunk as "current" accumulator
2. For each subsequent chunk, check if it can be absorbed
3. Absorption conditions:
   - Current content below `CHUNK_MIN_SIZE_TARGET`
   - Merging wouldn't exceed `CHUNK_SIZE`
   - Both chunks from same source document
4. If possible, merge and continue
5. If not, finalize current and start fresh with next chunk

**Key Design Decisions**:
- Forward-only merging (preserves natural document flow)
- No cross-document merging (preserves clear boundaries)
- Respects maximum size (never discards content)
- Metadata inheritance (merged chunks inherit from first chunk)
- Uses `\n\n` separator (preserves visual structure)

#### RAG Embedding Support

Change RAG embedding model in Admin Panel > Settings > Documents menu. Supports:
- Ollama models
- OpenAI models

#### Citations in RAG

The RAG feature includes citations for reference points, ensuring transparency and accountability in use of external sources.

#### File Context vs Builtin Tools

**File Context Capability**:

| File Context   | Behavior                                                      |
| -------------- | ------------------------------------------------------------- |
| Enabled        | Attached files processed via RAG. Content retrieved and injected.   |
| Disabled       | File processing completely skipped. No content extraction or injection. |

**Builtin Tools Capability**:

| Builtin Tools  | Behavior                                                                  |
| -------------- | ------------------------------------------------------------------------- |
| Enabled        | Model receives tools like `query_knowledge_bases`, `view_knowledge_file`, etc. |
| Disabled       | No builtin tools. Model works only with pre-injected context.             |

**Combination Matrix**:

| File Context   | Builtin Tools | Result                                                          |
| -------------- | ------------- | --------------------------------------------------------------- |
| Enabled        | Enabled       | Full Agentic Mode: RAG injected + autonomous querying            |
| Enabled        | Disabled      | Traditional RAG: Content injected, no autonomous tools           |
| Disabled       | Enabled       | Tools-Only Mode: No pre-injection, on-demand retrieval           |
| Disabled       | Disabled      | No File Processing: Files ignored                               |

#### Enhanced RAG Pipeline

Hybrid search sub-feature with:
- **BM25** keyword search
- **CrossEncoder** re-ranking
- Configurable relevance score thresholds

#### KV Cache Optimization

**Problem**: By default, RAG context is injected into user message, which shifts position during conversation and invalidates KV cache for each response.

**Solution**: Enable `RAG_SYSTEM_CONTEXT=True` environment variable
- Injects RAG context into system message instead
- System message stays at beginning, never changes position
- Provider can effectively cache processed context
- Follow-up questions benefit from instant responses and cost savings

**Recommended for**: Ollama, llama.cpp, OpenAI, Vertex AI users who frequently "chat with documents"

#### YouTube RAG Pipeline

Dedicated RAG pipeline for summarizing YouTube videos via video URLs. Incorporates video transcriptions directly into chats.

#### Document Parsing

Variety of parsers extract content from local and remote documents.

**Temporary Chat Limitations**: When using Temporary Chat, document processing is restricted to frontend-only operations to ensure privacy. Advanced backend parsing is disabled.

#### Google Drive Integration

Direct access to Drive files from chat interface. Requires:
1. Google Cloud project with Google Picker API and Google Drive API enabled
2. Environment variables: `GOOGLE_DRIVE_API_KEY` and `GOOGLE_DRIVE_CLIENT_ID`
3. Enable in Admin Panel > Settings > Documents > Google Drive

**Setup Steps**:
1. Create OAuth 2.0 client
2. Configure Authorized JavaScript origins & redirect URI
3. Note the Client ID
4. Enable Google Drive API and Google Picker API
5. Set app as Testing and add email to User List
6. Create API key with Website restrictions
7. Set up environment variables
8. Relaunch Open-WebUI instance

### Updating Open WebUI

Keeping Open WebUI updated ensures latest features, security patches, and bug fixes.

#### Before Updating

- **Backup your data** before major version updates
- **Check release notes** at `https://github.com/open-webui/open-webui/releases`
- **Clear browser cache** after updating
- **Multiple Workers**: If using `UVICORN_WORKERS > 1`:
  - Run updated container with `UVICORN_WORKERS=1` first, OR
  - Designate master worker with `ENABLE_DB_MIGRATIONS=True` on one instance

#### Manual Update

**Step 1: Stop and Remove Current Container**

```bash
docker rm -f open-webui
```

**Step 2: Pull Latest Docker Image**

```bash
docker pull ghcr.io/open-webui/open-webui:v0.8.3
```

**Step 3: Start Container with Updated Image**

```bash
docker run -d \
  -p 3000:8080 \
  -v open-webui:/app/backend/data \
  --name open-webui \
  --restart always \
  ghcr.io/open-webui/open-webui:v0.8.3
```

**Step 4: Verify Update**

```bash
# Check logs
docker logs open-webui

# Verify in browser
# 1. Navigate to http://localhost:3000
# 2. Clear browser cache
# 3. Hard refresh (Ctrl+F5)
# 4. Login and verify data
```

#### Using Image Tags

**For latest version**:
- `ghcr.io/open-webui/open-webui:v0.8.3`
- `ghcr.io/open-webui/open-webui:cuda`
- `ghcr.io/open-webui/open-webui:ollama`

**For production** (pinned versions):
- `ghcr.io/open-webui/open-webui:v0.8.3`
- `ghcr.io/open-webui/open-webui:v0.8.3-cuda`
- `ghcr.io/open-webui/open-webui:v0.8.3-ollama`

#### Persistent Login Sessions

To avoid being logged out after updates, set persistent `WEBUI_SECRET_KEY`:

```bash
docker run -d \
  -p 3000:8080 \
  -v open-webui:/app/backend/data \
  --name open-webui \
  --restart always \
  -e WEBUI_SECRET_KEY="your-secret-key-here" \
  ghcr.io/open-webui/open-webui:v0.8.3
```

**Generate secure key**:
```bash
openssl rand -hex 32
# or
python3 -c "import secrets; print(secrets.token_hex(32))"
```

#### Automated Updates

**Important**: Automated updates can break deployments if:
- New version has breaking changes
- Custom configurations become incompatible
- Database migrations fail during unattended updates

**Watchtower (Recommended Fork)**

The original `containrrr/watchtower` is not maintained. Use `nickfedor/watchtower` fork instead.

**One-time update**:
```bash
docker run --rm \
  --volume /var/run/docker.sock:/var/run/docker.sock \
  nickfedor/watchtower \
  --run-once open-webui
```

**Configuration options**:
| Variable                     | Description                      | Default      |
| ---------------------------- | -------------------------------- | ------------ |
| `WATCHTOWER_CLEANUP`         | Remove old images after update   | `false`      |
| `WATCHTOWER_INCLUDE_STOPPED` | Update stopped containers        | `false`      |
| `WATCHTOWER_SCHEDULE`        | Cron expression for schedule     | `0 0 0 * * *`|
| `WATCHTOWER_MONITOR_ONLY`    | Only notify, don't update        | `false`      |

**Monitor-only mode**:
```bash
docker run -d \
  --name watchtower \
  --volume /var/run/docker.sock:/var/run/docker.sock \
  -e WATCHTOWER_MONITOR_ONLY=true \
  -e WATCHTOWER_NOTIFICATIONS=email \
  -e WATCHTOWER_NOTIFICATION_EMAIL_TO=you@example.com \
  nickfedor/watchtower
```

**What's Up Docker (WUD)** - Web UI for visual monitoring and manual updates

**Diun** - Lightweight notification-only tool

#### Docker Volume Management

**Locate data**:
```bash
docker volume inspect open-webui
```

**Common locations**:
- Linux: `/var/lib/docker/volumes/open-webui/_data`
- Windows (WSL2): `\\wsl$\docker-desktop\mnt\docker-desktop-disk\data\docker\volumes\open-webui\_data`
- macOS: `~/Library/Containers/com.docker.docker/Data/vms/0/data/docker/volumes/open-webui/_data`

**Backup volume**:
```bash
docker run --rm \
  -v open-webui:/data \
  -v $(pwd):/backup \
  alpine tar czf /backup/openwebui-$(date +%Y%m%d_%H%M%S).tar.gz /data
```

**Restore volume**:
```bash
docker stop open-webui
docker run --rm \
  -v open-webui:/data \
  -v $(pwd):/backup \
  alpine sh -c "rm -rf /data/* && tar xzf /backup/openwebui-20241201_120000.tar.gz -C /"
docker start open-webui
```

**Clean up old images**:
```bash
docker image prune
docker image prune -a
```

#### Post-Update Checklist

- Open WebUI starts without errors
- Can access web interface
- Can login with existing credentials
- Chat history is intact
- Models configured correctly
- Custom settings preserved
- No JavaScript console errors
- Clear browser cache if needed

---

## API Complete Reference

### Authentication

API requests require authentication using the Bearer Token method. Retrieve your API key from **Settings > Account** in the Open WebUI, or utilize a JWT (JSON Web Token).

### Notable API Endpoints

### Retrieve All Models

- **Endpoint**: `GET /api/models`
- **Description**: Fetches all models created or added through Open WebUI.

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" http://localhost:3000/api/models
```

### Chat Completions

- **Endpoint**: `POST /api/chat/completions`
- **Description**: OpenAI API-compatible chat completion endpoint for models in Open WebUI, including Ollama, OpenAI, and Function models.

**Example:**
```bash
curl -X POST http://localhost:3000/api/chat/completions \
-H "Authorization: Bearer YOUR_API_KEY" \
-H "Content-Type: application/json" \
-d '{
"model": "llama3.1",
"messages": [{"role": "user", "content": "Why is the sky blue?"}]
}'
```

**Python Example:**
```python
import requests

def chat_with_model(token):
    url = 'http://localhost:3000/api/chat/completions'
    headers = {
        'Authorization': f'Bearer {token}',
        'Content-Type': 'application/json'
    }
    data = {
      "model": "granite3.1-dense:8b",
      "messages": [
        {
          "role": "user",
          "content": "Why is the sky blue?"
        }
      ]
    }
    response = requests.post(url, headers=headers, json=data)
    return response.json()
```

### Ollama API Proxy Support

If you want to interact directly with Ollama models—including for embedding generation or raw prompt streaming—Open WebUI offers a transparent passthrough to the native Ollama API via a proxy route.

- **Base URL**: `/ollama/<api>`

**Generate Completion (Streaming):**
```bash
curl http://localhost:3000/ollama/api/generate \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
  "model": "llama3.2",
  "prompt": "Why is the sky blue?"
}'
```

**List Available Models:**
```bash
curl http://localhost:3000/ollama/api/tags \
  -H "Authorization: Bearer YOUR_API_KEY"
```

**Generate Embeddings:**
```bash
curl -X POST http://localhost:3000/ollama/api/embed \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
  "model": "llama3.2",
  "input": ["Open WebUI is great!", "Let's generate embeddings."]
}'
```

### Retrieval Augmented Generation (RAG) API

### Uploading Files

Upload files to be used in RAG responses; the content is extracted and stored in a vector database.

- **Endpoint**: `POST /api/v1/files/`

**Curl Example:**
```bash
curl -X POST -H "Authorization: Bearer YOUR_API_KEY" -H "Accept: application/json" \
-F "file=@/path/to/your/file" http://localhost:3000/api/v1/files/
```

**Python Example:**
```python
import requests

def upload_file(token, file_path):
    url = 'http://localhost:3000/api/v1/files/'
    headers = {'Authorization': f'Bearer {token}', 'Accept': 'application/json'}
    files = {'file': open(file_path, 'rb')}
    response = requests.post(url, headers=headers, files=files)
    return response.json()
```

### Checking File Processing Status

Before adding a file to a knowledge base, verify that processing has completed.

- **Endpoint**: `GET /api/v1/files/{id}/process/status`

**Status Values:**

| Status     | Description                                                     |
| ---------- | --------------------------------------------------------------- |
| `pending`  | File is still being processed                                   |
| `completed`| Processing finished successfully                               |
| `failed`   | Processing failed (check `error` field for details)            |

**Python Example (Polling):**
```python
import requests
import time

def wait_for_file_processing(token, file_id, timeout=300, poll_interval=2):
    """Wait for a file to finish processing."""
    url = f'http://localhost:3000/api/v1/files/{file_id}/process/status'
    headers = {'Authorization': f'Bearer {token}'}

    start_time = time.time()
    while time.time() - start_time < timeout:
        response = requests.get(url, headers=headers)
        result = response.json()
        status = result.get('status')

        if status == 'completed':
            return result
        elif status == 'failed':
            raise Exception(f"File processing failed: {result.get('error')}")

        time.sleep(poll_interval)

    raise TimeoutError(f"File processing did not complete within {timeout} seconds")
```

### Adding Files to Knowledge Collections

Group uploaded files into a knowledge collection for reference in chats.

- **Endpoint**: `POST /api/v1/knowledge/{id}/file/add`

**Curl Example:**
```bash
curl -X POST http://localhost:3000/api/v1/knowledge/{knowledge_id}/file/add \
-H "Authorization: Bearer YOUR_API_KEY" \
-H "Content-Type: application/json" \
-d '{"file_id": "your-file-id-here"}'
```

**Python Example:**
```python
import requests

def add_file_to_knowledge(token, knowledge_id, file_id):
    url = f'http://localhost:3000/api/v1/knowledge/{knowledge_id}/file/add'
    headers = {
        'Authorization': f'Bearer {token}',
        'Content-Type': 'application/json'
    }
    data = {'file_id': file_id}
    response = requests.post(url, headers=headers, json=data)
    return response.json()
```

### Using Files and Collections in Chat Completions

**Using an Individual File:**
```bash
curl -X POST http://localhost:3000/api/chat/completions \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
        "model": "gpt-4-turbo",
        "messages": [
          {"role": "user", "content": "Explain the concepts in this document."}
        ],
        "files": [
          {"type": "file", "id": "your-file-id-here"}
        ]
      }'
```

**Using a Knowledge Collection:**
```bash
curl -X POST http://localhost:3000/api/chat/completions \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
        "model": "gpt-4-turbo",
        "messages": [
          {"role": "user", "content": "Provide insights on the historical perspectives covered in the collection."}
        ],
        "files": [
          {"type": "collection", "id": "your-collection-id-here"}
        ]
      }'
```

### Complete API Endpoint Tree

### Primary UI & Management Endpoints (`/api`)

**Chat Management (`/api/chats`)**
- `POST /api/chats/new` - Create new chat session
- `GET /api/chats` - List user's chat sessions
- `GET /api/chats/{chat_id}` - Retrieve specific chat with full history
- `POST /api/chats/{chat_id}` - Update entire chat object including messages
- `DELETE /api/chats/{chat_id}` - Delete chat session

**Chat Processing**
- `POST /api/chat/completions` - Unified chat completion endpoint (OpenAI-compatible)
- `POST /api/chat/completed` - Completion finalization for post-processing

**Model Management (`/api/models`)**
- `GET /api/models` - List available models and metadata
- `POST /api/models/pull` - Pull or download models from registries

**Prompt Management (`/api/prompts`)**
- `GET /api/prompts` - List prompt templates
- `POST /api/prompts` - Create new prompt templates

**Administrative Functions**
- User Management (`/api/users`) - List, update, and delete users (admin only)
- System Management (`/api/system`) - Reload configurations and status checks (admin only)

### File & RAG Management (`/api/v1`)

**File Operations (`/api/v1/files`)**
- `POST /api/v1/files` - Upload files via multipart form data
- `GET /api/v1/files` - List uploaded files
- `GET /api/v1/files/{id}/content` - Download raw file content
- `DELETE /api/v1/files/{id}` - Delete uploaded files

**Knowledge Base Operations (`/api/v1/knowledge`)**
- `POST /api/v1/knowledge/create` - Create new knowledge bases
- `GET /api/v1/knowledge` - List available knowledge bases
- `GET /api/v1/knowledge/{id}` - Get knowledge base details
- `POST /api/v1/knowledge/{id}/file/add` - Add files to knowledge bases
- `DELETE /api/v1/knowledge/{id}/delete` - Delete knowledge bases

### Integration Endpoints

**Ollama Passthrough (`/ollama`)**
- `POST /ollama/api/generate` - Native Ollama generate endpoint
- `POST /ollama/api/chat` - Native Ollama chat endpoint
- `GET /ollama/api/tags` - List available Ollama models
- `POST /ollama/api/embed` - Generate embeddings using Ollama

**OpenAI Compatibility (`/v1`)**
- `GET /v1/models` - OpenAI-style model listing
- `POST /v1/chat/completions` - OpenAI-compatible chat completions

### Tools and Functions Endpoints

**Tools (`/api/v1/tools`)**
- `GET /api/v1/tools/` - List available tools
- `GET /api/v1/tools/id/{id}` - Get specific tool details

**Functions (`/api/v1/functions`)**
- `GET /api/v1/functions/` - List available functions

## Detailed API Endpoint Reference

### File Download Endpoint

#### Download File Content

```http
GET /api/v1/files/{id}/content
Authorization: Bearer <token>
```

**Response:**
- Content-Type: based on file type
- Body: Raw file content

**Errors:**
- `404 Not Found`: File doesn't exist
- `403 Forbidden`: No access to file

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/files/{file_id}/content \
  -o downloaded_file.pdf
```

### Chat Management Endpoints

#### Create New Chat

```http
POST /api/v1/chats/new
Authorization: Bearer <token>
Content-Type: application/json
```

**Request Body:**
```json
{
  "chat": {
    "title": "New Chat",
    "models": ["model-id"],
    "messages": [],
    "history": {}
  }
}
```

**Response (201 Created):**
Returns the created chat object with generated ID.

**Example:**
```bash
curl -X POST http://localhost:3000/api/v1/chats/new \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "chat": {
      "title": "Research Discussion",
      "models": ["llama3.2"],
      "messages": [],
      "history": {}
    }
  }'
```

#### List User Chats

```http
GET /api/v1/chats
Authorization: Bearer <token>
```

**Response (200 OK):**
```json
{
  "chats": [
    {
      "id": "chat-id",
      "title": "Chat Title",
      "created_at": "2024-01-15T10:30:00Z",
      "updated_at": "2024-01-15T10:35:00Z"
    }
  ]
}
```

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/chats | jq .
```

#### Get Chat by ID

```http
GET /api/v1/chats/{id}
Authorization: Bearer <token>
```

**Response:** Returns the complete chat object including all messages and metadata.

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/chats/{chat_id} | jq .
```

#### Update Chat

```http
POST /api/v1/chats/{id}
Authorization: Bearer <token>
Content-Type: application/json
```

**Request Body:** Complete chat object with updates

**Example:**
```bash
curl -X POST http://localhost:3000/api/v1/chats/{chat_id} \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "title": "Updated Chat Title",
    "messages": [...],
    "history": {...}
  }'
```

#### Delete Chat

```http
DELETE /api/v1/chats/{id}
Authorization: Bearer <token>
```

**Response:** `204 No Content`

**Example:**
```bash
curl -X DELETE -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/chats/{chat_id}
```

### Tools and Functions Endpoints

#### List Tools

```http
GET /api/v1/tools/
Authorization: Bearer <token>
```

**Response (200 OK):**
```json
{
  "tools": [
    {
      "id": "tool-id",
      "name": "Tool Name",
      "description": "Tool description",
      "enabled": true
    }
  ]
}
```

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/tools/ | jq .
```

#### Get Tool Details

```http
GET /api/v1/tools/id/{id}
Authorization: Bearer <token>
```

**Response:** Returns detailed information about a specific tool including its schema and configuration.

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/tools/id/{tool_id} | jq .
```

#### List Functions

```http
GET /api/v1/functions/
Authorization: Bearer <token>
```

**Response:** Similar structure to tools endpoint, listing available functions.

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/api/v1/functions/ | jq .
```

### Ollama Proxy Endpoints

#### Ollama Generate

```http
POST /ollama/api/generate
Authorization: Bearer <token>
Content-Type: application/json
```

**Request Body:**
```json
{
  "model": "llama3.2",
  "prompt": "Why is the sky blue?",
  "stream": false
}
```

**Example:**
```bash
curl -X POST http://localhost:3000/ollama/api/generate \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.2",
    "prompt": "Why is the sky blue?",
    "stream": false
  }'
```

#### Ollama Chat

```http
POST /ollama/api/chat
Authorization: Bearer <token>
Content-Type: application/json
```

**Request Body:**
```json
{
  "model": "llama3.2",
  "messages": [
    {"role": "user", "content": "Hello!"}
  ],
  "stream": false
}
```

**Example:**
```bash
curl -X POST http://localhost:3000/ollama/api/chat \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.2",
    "messages": [
      {"role": "user", "content": "Hello!"}
    ],
    "stream": false
  }'
```

#### List Ollama Models

```http
GET /ollama/api/tags
Authorization: Bearer <token>
```

**Example:**
```bash
curl -H "Authorization: Bearer YOUR_API_KEY" \
  http://localhost:3000/ollama/api/tags | jq .
```

#### Generate Embeddings (Ollama)

```http
POST /ollama/api/embed
Authorization: Bearer <token>
Content-Type: application/json
```

**Request Body:**
```json
{
  "model": "llama3.2",
  "input": ["Text to embed", "Another text"]
}
```

**Example:**
```bash
curl -X POST http://localhost:3000/ollama/api/embed \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "llama3.2",
    "input": ["Text to embed", "Another text"]
  }'
```

## Error Handling

### HTTP Status Codes

| Code | Meaning | Common Causes |
|------|---------|---------------|
| `200 OK` | Success | Request successful |
| `201 Created` | Resource created | New chat/file created |
| `204 No Content` | Success (no body) | Deletion successful |
| `400 Bad Request` | Invalid request | Malformed JSON, missing fields |
| `401 Unauthorized` | Not authenticated | Missing or invalid token |
| `403 Forbidden` | No permission | Valid token but insufficient rights |
| `404 Not Found` | Resource missing | Invalid ID, deleted resource |
| `422 Unprocessable Entity` | Validation error | Schema validation failed |
| `429 Too Many Requests` | Rate limited | Too many requests |
| `500 Internal Server Error` | Server error | Upstream or internal failure |

### Error Response Format

```json
{
  "detail": "Human-readable error message"
}
```

### Common Error Scenarios

**File Upload Errors:**
- `400`: File too large (check RAG_FILE_MAX_SIZE)
- `400`: Unsupported file type (check RAG_ALLOWED_FILE_EXTENSIONS)
- `408`: Processing timeout

**RAG Errors:**
- `500`: Embedding model not available
- `500`: Vector store connection failed
- `422`: Document parsing failed

**Chat Completion Errors:**
- `404`: Model not found
- `400`: Context length exceeded
- `429`: Rate limit exceeded

## Swagger Documentation Links

Access API documentation for Open WebUI services:

| Application | Documentation Path |
|-------------|-------------------|
| Main | `/docs` |
| WebUI | `/api/v1/docs` |
| Ollama | `/ollama/docs` |
| OpenAI | `/openai/docs` |
| Images | `/images/api/v1/docs` |
| Audio | `/audio/api/v1/docs` |
| RAG | `/retrieval/api/v1/docs` |

---

## Complete Troubleshooting Guide

## General Troubleshooting Tips

Encountering issues? Start with these important steps:

- Make sure you're using the **latest version** of the software
- **Check for Configuration Persistence:** Open WebUI prioritizes settings stored in its internal database over environment variables for certain settings (marked as `PersistentConfig`)

With this project constantly evolving, updates and fixes are regularly added. Keeping your software up-to-date is crucial.

**Community Support:** If you still face problems after updating, join our vibrant community on Discord.

## Web Search Troubleshooting

### Web Search Fails Behind HTTP Proxy

If you're running Open WebUI behind an HTTP proxy, you might notice that web search queries succeed, but the subsequent content fetching fails with errors like:
- `[Errno -3] Temporary failure in name resolution`
- `Connection timeout to host`
- `The content provided is empty`

**Solution:**

Navigate to: Admin Panel > Settings > Web Search
- Enable "Trust Proxy Environment"
- Save changes

Alternatively, set the environment variable:
```bash
WEB_SEARCH_TRUST_ENV=True
```

### 403 Forbidden Errors from SearXNG

If you're using SearXNG and seeing `403 Client Error: Forbidden` in your logs, the JSON format is not enabled.

**Solution:**

Edit your SearXNG `settings.yml` and add `json` to the formats list:
```yaml
search:
  formats:
    - html
    - json
```

Restart SearXNG after making this change.

### Empty Content or Poor Results

**Solutions:**

- **Increase context length:** Web pages often contain 4,000-8,000+ tokens. Increase to 16384+ tokens in Admin Panel > Models > Settings > Advanced Parameters
- **Check result count:** Adjust `WEB_SEARCH_RESULT_COUNT`
- **Try different loaders:** Configure `WEB_LOADER_ENGINE` to use `playwright` for JavaScript-heavy sites or `firecrawl/tavily` for better extraction

### Environment Variables Reference

| Variable | Description |
|----------|-------------|
| `WEB_SEARCH_TRUST_ENV` | Enable proxy support for content fetching |
| `WEB_SEARCH_RESULT_COUNT` | Number of search results to fetch |
| `WEB_SEARCH_CONCURRENT_REQUESTS` | Concurrent requests to search engine |
| `WEB_LOADER_CONCURRENT_REQUESTS` | Concurrent page fetches |
| `WEB_LOADER_ENGINE` | Content extraction engine |

## Audio Troubleshooting

### Microphone Access Issues

For security reasons, accessing the microphone is restricted to pages served over HTTPS or locally from localhost.

**Solutions for Non-HTTPS Connections:**

1. **Set Up HTTPS (Recommended):** Configure your server with a reverse proxy like Nginx or Caddy with Let's Encrypt certificates

2. **Temporary Browser Flags (Use with caution):**

**Chromium-based Browsers:**
- Open `chrome://flags/#unsafely-treat-insecure-origin-as-secure`
- Enter your non-HTTPS address
- Restart the browser

**Firefox:**
- Open `about:config`
- Modify `dom.securecontext.allowlist`
- Add your IP addresses separated by commas

### Text-to-Speech (TTS) Issues

#### TTS Loading Forever / Not Working

**Symptoms:**
- TTS keeps loading forever
- Container logs show: `RuntimeError: Dataset scripts are no longer supported`

**Solutions:**
- Temporary fix: `docker exec open-webui bash -lc "pip install datasets==3.6.0" && docker restart open-webui`
- Permanent fix: Add `EXTRA_PIP_PACKAGES=datasets==3.6.0` to docker-compose.yml

### Speech-to-Text (STT) Issues

#### Whisper STT Not Working / Compute Type Error

**Symptoms:**
- Error message: `Error transcribing chunk: Requested int8 compute type, but the target device or backend do not support efficient int8 computation`

**Solutions:**

1. **Upgrade to the Latest Version (Recommended)**

2. **Manually Set Compute Type:**
```bash
WHISPER_COMPUTE_TYPE=float16
```

3. **Switch to the Standard Image:**
```bash
# Instead of: ghcr.io/open-webui/open-webui:cuda
# Use: ghcr.io/open-webui/open-webui:v0.8.3
```

**Available Compute Types:**

| Compute Type | Best For | Notes |
|--------------|----------|-------|
| `int8` | CPU (default) | Fastest, but doesn't work on older GPUs |
| `float16` | CUDA/GPU (recommended) | Best balance of speed and compatibility |
| `int8_float16` | GPU with hybrid precision | Uses int8 for weights, float16 for computation |
| `float32` | Maximum compatibility | Slowest, but works on all hardware |

## RAG Troubleshooting

### Common RAG Issues

#### 1. The Model "Can't See" Your Content

This is the most common problem—and it's typically caused by issues during your content ingestion process.

**Solution:** Check your content extraction settings

Navigate to: Admin Settings > Documents.
Make sure you're using a robust content extraction engine such as:
- Apache Tika
- Docling
- Custom extractors (depending on your document types)

#### 2. Only a Small Part of the Document is Being Used

Open WebUI aggressively trims down the retrieved content to fit within the assumed available space.

**Solutions:**

- **Enable "Bypass Embedding and Retrieval"** — This sends full content directly without applying strict retrieval filters
- **Toggle on "Full Context Mode"** — This injects more comprehensive content into the model prompt

#### 3. Token Limit is Too Short

By default, many models are limited to a 2048-token context window.

**Solutions:**

- **For Ollama Models:** Extend the model's context length in Admin Panel > Models > Settings > Advanced Parameters (increase to 8192+ or ideally beyond 16000 tokens)
- **For OpenAI and Other Integrated Models:** Ensure you're using a model with sufficient context length
- **Use an external LLM** with larger context capacity (GPT-4, GPT-4o, Claude 3, Gemini 1.5, or Mixtral with 8k+ context)

#### 4. Embedding Model is Low-Quality or Mismatched

**Solution:**

Change to a high-quality embedding model (e.g., all-MiniLM-L6-v2, Instructor X, or OpenAI embeddings)
- Go to: Admin Settings > Documents > Embedding Model
- Save the embedding model again—even if it's already selected
- After changing the model, **Reindex all existing documents**

#### 5. 400: 'NoneType' object has no attribute 'encode'

**Cause:** Your embedding model isn't set up properly.

**Solution:**
- Go to: Admin Settings > Documents > Embedding Model
- Save the embedding model again
- If using a remote/external embedding tool, make sure it's running and accessible

#### 6. Upload Limits and Restrictions

**Chat Uploads:**
- Max File Size: `RAG_FILE_MAX_SIZE` (default: Unlimited)
- Max File Count: `RAG_FILE_MAX_COUNT` (default: Unlimited)
- Allowed File Extensions: `RAG_ALLOWED_FILE_EXTENSIONS` (default: All)

**Folder Uploads:**
- Subject to `FOLDER_MAX_FILE_COUNT` (defaults to 100)

**Knowledge Base Uploads:**
- File Limit: Subject to `RAG_FILE_MAX_SIZE`, but NOT `RAG_FILE_MAX_COUNT`
- RAG Enforcement: All files automatically indexed

#### 7. Fragmented or Tiny Chunks

**Solution:**
- Go to Admin Settings > Documents
- Increase the **Chunk Min Size Target**
- Setting to ~1000 (or 50-60% of your `CHUNK_SIZE`) will merge small fragments

#### 8. Slow Follow-up Responses (KV Cache Invalidation)

**Solution:**
```bash
RAG_SYSTEM_CONTEXT=True
```
This injects the RAG context into the system message, which stays at a fixed position at the start of the conversation.

#### 9. CUDA Out of Memory During Embedding

**Solutions:**

- Isolate Embedding to a Different GPU: Set `CUDA_VISIBLE_DEVICES`
- Reduce Embedding Batch Size: Lower `RAG_EMBEDDING_BATCH_SIZE` (e.g., from 32 to 8 or 4)
- Enable Expandable Segments: `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`
- Restart Between Large Ingestion Jobs
- Use Smaller Embedding Models

## OAuth/SSO Troubleshooting

### Common OAuth/SSO Issues

#### 1. WebUI URL Not Configured in Admin Panel

**Solution:**

Navigate to: Admin Settings > General
- Ensure your WebUI URL field is filled in and points to your deployed instance (e.g., `https://yourwebui.yourdomain.com`)

#### 2. Incorrect Environment Variable Configuration

**Incorrect Variables People Often Use:**
- `OIDC_CONFIG` - Use `OPENID_PROVIDER_URL` instead
- `WEBUI_OIDC_CLIENT_ID` - Use `OAUTH_CLIENT_ID` instead
- `WEBUI_ENABLE_SSO` - Use `ENABLE_OAUTH_SIGNUP` instead
- `WEBUI_AUTH_TYPE` - This doesn't exist
- `OPENID_CLIENT_ID` - Use `OAUTH_CLIENT_ID` instead
- `OPENID_CLIENT_SECRET` - Use `OAUTH_CLIENT_SECRET` instead

**Correct OIDC Variables:**
```bash
# Required for OIDC
OAUTH_CLIENT_ID=your_client_id
OAUTH_CLIENT_SECRET=your_client_secret
OPENID_PROVIDER_URL=https://your-provider/.well-known/openid-configuration
ENABLE_OAUTH_SIGNUP=true

# Optional but recommended
OAUTH_PROVIDER_NAME=Your Provider Name
OAUTH_SCOPES=openid email profile
OPENID_REDIRECT_URI=https://your-domain/oauth/oidc/callback
```

**Correct Microsoft Variables:**
```bash
MICROSOFT_CLIENT_ID=your_client_id
MICROSOFT_CLIENT_SECRET=your_client_secret
MICROSOFT_CLIENT_TENANT_ID=your_tenant_id
OPENID_PROVIDER_URL=https://login.microsoftonline.com/YOUR_TENANT_ID/v2.0/.well-known/openid-configuration
ENABLE_OAUTH_SIGNUP=true
```

#### 3. Missing Required Variables

`OPENID_PROVIDER_URL` is mandatory for OIDC

- **Microsoft Entra ID:** `https://login.microsoftonline.com/YOUR_TENANT_ID/v2.0/.well-known/openid-configuration`
- **Google:** `https://accounts.google.com/.well-known/openid-configuration`
- **Authentik:** `https://your-authentik-domain/application/o/your-app-name/.well-known/openid-configuration`

#### 4. Persistent Configuration Conflicts

New Issue: OAuth settings are stored in the database after the first launch when `ENABLE_OAUTH_PERSISTENT_CONFIG=true` (default).

**Solutions:**

- **For Development/Testing:** Set `ENABLE_OAUTH_PERSISTENT_CONFIG=false`
- **For Production:** Either configure settings through Admin Panel OR temporarily set `ENABLE_OAUTH_PERSISTENT_CONFIG=false`, restart to apply new env vars, then set back to true
- **Fresh Start:** Delete the database volume and restart

#### 5. Server-Side Caching (A Hidden Trouble Spot!)

**Solutions:**

In your NGINX configuration, exclude these endpoints from server-side caching:
- `/api`, `/oauth`, `/callback`, `/login`, `/ws`, `/websocket`

**Example NGINX Configuration:**
```nginx
location ~* ^/(api|oauth|callback|login|ws|websocket) {
    proxy_no_cache 1;
    proxy_cache_bypass 1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
    proxy_http_version 1.1;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";
    proxy_read_timeout 3600s;
    proxy_send_timeout 3600s;
    proxy_set_header Accept-Encoding "";
}
```

## Multi-Replica / High Availability / Concurrency

### Core Requirements Checklist

Before troubleshooting, ensure your deployment meets these absolute requirements:

- **Shared Secret Key:** `WEBUI_SECRET_KEY` MUST be identical on all replicas
- **External Database:** You MUST use an external PostgreSQL database (see `DATABASE_URL`). SQLite is NOT supported
- **Redis for WebSockets:** `ENABLE_WEBSOCKET_SUPPORT=True` and `WEBSOCKET_MANAGER=redis` with valid `WEBSOCKET_REDIS_URL`
- **Shared Storage:** A persistent volume (RWX) is critical for RAG and generated images
- **External Vector Database (Recommended):** Use dedicated Vector DB (e.g., PGVector, Milvus, Qdrant)
- **Database Session Sharing (Optional):** `DATABASE_ENABLE_SESSION_SHARING=True` for high concurrency

### Common Issues

#### 1. Login Loops / 401 Unauthorized Errors

**Cause:** Each replica is using a different `WEBUI_SECRET_KEY`

**Solution:**
```yaml
env:
  - name: WEBUI_SECRET_KEY
    value: "your-super-secure-static-key-here"
```

#### 2. WebSocket 403 Errors / Connection Failures

**Solutions:**

**Configure CORS:**
```bash
CORS_ALLOW_ORIGIN="https://chat.yourdomain.com;http://chat.yourdomain.com;https://yourhostname;http://localhost:3000"
```

**Enable Redis for WebSockets:**
```bash
ENABLE_WEBSOCKET_SUPPORT=True
WEBSOCKET_MANAGER=redis
WEBSOCKET_REDIS_URL=redis://your-redis-host:6379/0
```

#### 3. Database Corruption / "Locked" Errors

**Cause:** Using SQLite with multiple replicas

**Solution:** Migrate to PostgreSQL
```bash
DATABASE_URL=postgresql://user:password@postgres-host:5432/openwebui
```

## Server Connectivity Issues

### HTTPS, TLS, CORS & WebSocket Issues

#### Common Symptoms

- Empty responses like "{}" in the chat
- Errors like "Unexpected token 'd', "data: {"id"... is not valid JSON"
- Garbled markdown during streaming
- WebSocket connection failures
- Login problems or session issues
- CORS errors in browser developer tools
- Mixed content warnings

#### Required Configuration for HTTPS & Reverse Proxies

```bash
# Set this to your actual domain BEFORE FIRST STARTUP
WEBUI_URL=https://your-open-webui-domain.com

# CORS configuration - CRITICAL for WebSocket functionality
CORS_ALLOW_ORIGIN="https://yourdomain.com;http://yourdomain.com;https://yourip;http://localhost:3000"

# Cookie security settings for HTTPS
WEBUI_SESSION_COOKIE_SECURE=true
WEBUI_AUTH_COOKIE_SECURE=true

# For OAuth/SSO
WEBUI_SESSION_COOKIE_SAME_SITE=lax
WEBUI_AUTH_COOKIE_SAME_SITE=lax

# WebSocket support (if using Redis)
ENABLE_WEBSOCKET_SUPPORT=true
WEBSOCKET_MANAGER=redis
WEBSOCKET_REDIS_URL=redis://redis:6379/1
```

### Garbled Markdown / Streaming Response Corruption

#### Cause: Nginx Proxy Buffering

**Solution: Disable Proxy Buffering**

```nginx
location / {
    proxy_pass http://your-open-webui-upstream;

    # CRITICAL: Disable buffering for SSE streaming
    proxy_buffering off;
    proxy_cache off;

    # WebSocket support
    proxy_http_version 1.1;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";

    # Standard proxy headers
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
}
```

### Model List Loading Issues

**Solution 1: Adjust the Timeout**
```bash
AIOHTTP_CLIENT_TIMEOUT_MODEL_LIST=3
```

**Solution 2: Fix or Remove Unreachable Endpoints**

**Solution 3: Recover from Database-Persisted Bad Configuration**
```bash
# Option A: Reset configuration on startup
RESET_CONFIG_ON_START=true

# Option B: Always use environment variables
ENABLE_PERSISTENT_CONFIG=false
```

## Image Generation Troubleshooting

### General Issues

**Image Not Generating:**
1. Check Images settings in Admin Panel > Settings > Images
2. Verify API Key and Base URL
3. Ensure selected model is available

### ComfyUI Issues

**Incompatible Workflow / JSON Errors:**
- In ComfyUI: Click "Settings" (gear icon)
- Enable "Enable Dev mode Options"
- Click "Save (API Format)" in the menu

**Image Editing / Image Variation Fails:**
- Custom workflow must have nodes configured to accept an input image
- Check default "Image Editing" workflow for required node structure

### Automatic1111 Issues

**Connection Refused / "Api Not Found":**
- Ensure Automatic1111 is running with `--api` flag

**Docker Connectivity:**
```bash
# Use host.docker.internal
http://host.docker.internal:7860
```

## Browser Compatibility

### Supported Browsers

Open WebUI's core functionality specifically depends on these browser versions:

- **Chrome 111** (Released March 2023)
- **Safari 16.4** (Released March 2023)
- **Firefox 128** (Released July 2024)

Note: If you experience any issues, ensure your browser is up to date or try an alternative supported browser.

## Reset Admin Password

### For Docker Deployments

**Step 1: Generate a New Password Hash**
```bash
htpasswd -bnBC 10 "" your-new-password | tr -d ':\n'
```

**Step 2: Update the Password in Docker**
```bash
docker run --rm -v open-webui:/data alpine/socat EXEC:"bash -c 'apk add sqlite && echo UPDATE auth SET password='\''HASH'\'' WHERE email='\''admin@example.com'\''; | sqlite3 /data/webui.db'", STDIO
```

### For Local Installations

**Step 1: Generate a New Password Hash**
```bash
htpasswd -bnBC 10 "" your-new-password | tr -d ':\n'
```

**Step 2: Update the Password Locally**
```bash
sqlite3 backend/data/webui.db "UPDATE auth SET password='HASH' WHERE email='admin@example.com';"
```

## Optimization, Performance & RAM Usage

### Performance Tuning (Speed & Responsiveness)

#### 1. Dedicated Task Models

Use a very fast, small, and cheap NON-REASONING model for background tasks (title generation, tagging, autocomplete).

**Configuration:**
- Admin Panel > Settings > Interface
- **Task Model (External):** Used when chatting with external model
- **Task Model (Local):** Used when chatting with local model

**Best Options:**
- External/Cloud: gpt-5-nano, gemini-2.5-flash-lite, llama-3.1-8b-instant
- Local: qwen3:1b, gemma3:1b, llama3.2:3b

#### 2. Caching & Latency Optimization

**Model Caching:**
```bash
ENABLE_BASE_MODELS_CACHE=True
MODELS_CACHE_TTL=300
```

**Search Query Caching:**
```bash
ENABLE_QUERIES_CACHE=True
```

**KV Cache Optimization (RAG Performance):**
```bash
RAG_SYSTEM_CONTEXT=True
```

### Database Optimization

#### PostgreSQL (Mandatory for Scale)

```bash
DATABASE_URL=postgres://user:password@localhost:5432/webui
```

#### Chat Saving Strategy

```bash
ENABLE_REALTIME_CHAT_SAVE=False  # Default - DO NOT ENABLE in production
```

#### Database Session Sharing

```bash
DATABASE_ENABLE_SESSION_SHARING
```

- **SQLite:** Keep disabled (default)
- **PostgreSQL with adequate resources:** Consider enabling for improved performance

#### Connection Pool Sizing

```bash
DATABASE_POOL_SIZE=15
DATABASE_POOL_MAX_OVERFLOW=20
```

### Resource Efficiency (Reducing RAM)

#### 1. Offload Auxiliary Models

**RAG Embeddings:**
```bash
RAG_EMBEDDING_ENGINE=openai  # Offload completely
```

**Speech-to-Text (STT):**
```bash
AUDIO_STT_ENGINE=webapi  # Zero server load
```

#### 2. Disable Unused Features

```bash
ENABLE_IMAGE_GENERATION=False
ENABLE_CODE_INTERPRETER=False
```

#### 3. Disable Background Tasks

- Autocomplete: `ENABLE_AUTOCOMPLETE_GENERATION=False` (High Impact)
- Follow-up Questions: `ENABLE_FOLLOW_UP_GENERATION=False`
- Title Generation: `ENABLE_TITLE_GENERATION=False`
- Tag Generation: `ENABLE_TAGS_GENERATION=False`

---

## Tutorials

### Tips & Tricks

### RAG Tutorial: Configuring RAG with Open WebUI Documentation

This tutorial demonstrates how to customize Open WebUI for your specific use case.

#### What is RAG?

Retrieval-Augmented Generation (RAG) combines LLMs with retrieved knowledge from external sources. The system retrieves relevant data from uploaded documents or knowledge bases, enhancing the quality and accuracy of responses.

#### Step-by-Step Setup: Open WebUI Documentation as Knowledge Base

1. **Download the Documentation:** `https://github.com/open-webui/docs/archive/refs/heads/main.zip`
2. **Extract the Files**
3. **Locate the Markdown Files** (*.md and *.mdx extensions)
4. **Create a Knowledge Base:** Navigate to Workspace > Knowledge > + Create a Knowledge Base
5. **Upload the Files**

#### Create and Configure the Model

1. Navigate to Workspace > Models > + Add New Model
2. Configure the Model:
   - Name: Open WebUI
   - Base Model: Select appropriate model
   - Knowledge Source: Select Open WebUI Documentation
3. Save the Model

#### Examples and Usage

Query the Open WebUI Documentation Model with questions like:
- "How do I configure environment variables?"
- "How do I update Open WebUI using Docker?"

### Dual OAuth Configuration (Microsoft & Google)

This tutorial covers configuring both Microsoft and Google OAuth simultaneously.

### SQLite Database Overview

Understanding how Open WebUI uses SQLite for data storage and management.

### Optimization, Performance & RAM Usage

Comprehensive guide to optimizing Open WebUI performance for different deployment scenarios.

#### Recommended Configuration Profiles

**Profile 1: Maximum Privacy (Weak Hardware/RPi)**
- Target: 100% Local, Raspberry Pi / <4GB RAM
- Embeddings: Default (SentenceTransformers)
- Audio: `AUDIO_STT_ENGINE=webapi`
- Task Model: Disable or use tiny model
- Scaling: Keep default `THREAD_POOL_SIZE` (40)
- Disable: Image Gen, Code Interpreter, Autocomplete, Follow-ups
- Database: SQLite is fine

**Profile 2: Single User Enthusiast**
- Target: Max Quality & Speed, Local + External APIs
- Embeddings: `RAG_EMBEDDING_ENGINE=openai`
- Task Model: gpt-5-nano or llama-3.1-8b-instant

## Integrations

### Azure OpenAI with EntraID

Complete guide for integrating Azure OpenAI services with Microsoft Entra ID authentication.

### Run DeepSeek R1 Dynamic 1.58-bit with Llama.cpp

Tutorial for running optimized DeepSeek R1 models using llama.cpp.

### Backend-Controlled, UI-Compatible API Flow

**Prerequisites:**
- Running Open WebUI instance
- API authentication token
- REST API knowledge

**7-Step Process:**

1. **Create chat with user message**
```bash
POST /api/chats/new
```

2. **Enrich chat response with assistant message**
```bash
POST /api/chats/{chat_id}
```

3. **Update chat with assistant message**
```bash
POST /api/chats/{chat_id}
```

4. **Trigger assistant completion**
```bash
POST /api/chat/completions
```

5. **Wait for response completion**

6. **Complete assistant message**
```bash
POST /api/chat/completed
```

7. **Fetch final chat**
```bash
GET /api/chats/{chat_id}
```

### Local LLM Setup with IPEX-LLM on Intel GPU

Guide for setting up Open WebUI with Intel GPU acceleration using IPEX-LLM.

### Continue.dev VS Code Extension with Open WebUI

**Download extension** from Visual Studio Marketplace

**Configure config.yaml:**
```yaml
name: Local Assistant
models:
  - name: LLama3
    provider: openai
    model: Meta-Llama-3-8B-Instruct-Q4_K_M.gguf
    apiBase: http://localhost:3000/api
    apiKey: YOUR_OPEN_WEBUI_API_KEY
    roles:
      - chat
      - edit
```

### Setting up with Custom CA Store

Guide for configuring Open WebUI with custom Certificate Authority stores.

### Browser Search Engine Integration

**Setup:**

1. Set `WEBUI_URL` environment variable
2. Add Open WebUI as custom search engine in Chrome/Firefox

**URL format:**
```
https://<your-open-webui-url>/?q=%s
```

**Optional specific model:**
```
https://<your-open-webui-url>/?models=<model_id>&q=%s
```

### Monitor your LLM requests with Helicone

Integration guide for monitoring Open WebUI API requests using Helicone.

### Monitoring and Debugging with Langfuse

Guide for integrating Langfuse observability platform.

### LibreTranslate Integration

How to integrate LibreTranslate for translation services.

### Redis Websocket Support

Configuration guide for using Redis to manage WebSocket connections in multi-instance deployments.

### Integrate with Amazon Bedrock

Complete guide for AWS Bedrock integration with Open WebUI.

### Integrate with MiniMax M2.1

Integration tutorial for MiniMax M2.1 API.

### Integrate with OneDrive & SharePoint

Guide for integrating Microsoft OneDrive and SharePoint for document access.

### Okta OIDC SSO Integration

Complete Okta Single Sign-On configuration guide.

### Notion MCP Integration

Integration guide for connecting Notion via Model Context Protocol.

### Jupyter Notebook Integration

Using Open WebUI within Jupyter Notebooks.

### Firefox AI Chatbot Sidebar

**Enable AI Chatbot in Firefox:**
- Enable in Firefox Labs or about:config
- Configure browser.ml.chat.provider with Open WebUI URL

**URL Parameters:**
- `models/model`: Specify model(s)
- `youtube`: Video transcription
- `web-search`: Enable web search
- `tools/tool-ids`: Activate specific tools
- `call`: Video call overlay
- `q`: Initial query
- `temporary-chat`: Temporary session

### Azure AD Domain Services (LDAPS) Integration

Guide for integrating with Azure Active Directory Domain Services via LDAPS.

### Iterm2 AI Integration

Using Open WebUI AI within iTerm2 terminal.

## Maintenance

### Exporting and Importing Database

Guide for database backup and migration.

### Switching to S3 Storage

Configuration guide for using AWS S3 or S3-compatible storage for file storage.

### Backups

Comprehensive backup strategy for Open WebUI deployments.

## HTTPS Setup

### HTTPS using Caddy

**Prerequisites:**
- Docker installed
- Domain name configured

**Create docker-compose.yml:**
```yaml
version: "3.8"

services:
  open-webui:
    image: ghcr.io/open-webui/open-webui:v0.8.3
    container_name: open-webui
    ports:
      - 3000:8080
    environment:
      - WEBUI_URL=https://your-domain.com
    volumes:
      - open-webui:/app/backend/data
    restart: always

  caddy:
    image: caddy:latest
    container_name: caddy
    ports:
      - 80:80
      - 443:443
    volumes:
      - ./Caddyfile:/etc/caddy/Caddyfile
      - caddy_data:/data
      - caddy_config:/config
    restart: always

volumes:
  open-webui:
  caddy_data:
  caddy_config:
```

**Configure Caddyfile:**
```text
your-domain.com {
    reverse_proxy open-webui:8080
}
```

**Test HTTPS at** https://your-domain.com

## Offline Mode

Guide for running Open WebUI completely offline without external API dependencies.

---


## About & Community

## Mission & Vision

Open WebUI is an open-source AI platform founded by **Tim J. Baek** with the mission to democratize AI through an accessible, user-friendly, and customizable interface.

**Core Values:**
- Open-source first
- Community-driven development
- Privacy-focused design
- Universal compatibility

## Team

**Founder:** Tim J. Baek

Open WebUI follows a community governance model with contributions from developers worldwide.

## Licensing

Open WebUI uses a **modified BSD-3 license** with branding protection (added in v0.6.6+):

**Permitted:**
- Self-hosting for any purpose
- Modification and customization
- Commercial use

**Restriction:**
- Enterprise entities must remove "Open WebUI" branding when using modified versions

## Contributing

### How to Contribute

1. **Testing Development Branch** - Run dev builds and report issues
2. **Code Contributions** - Submit pull requests following BSD-3 license
3. **Documentation** - Improve documentation and tutorials
4. **Translation** - Help translate the interface
5. **Accessibility** - Improve accessibility features

### Development Setup

See [Local Development Guide](https://docs.openwebui.com/getting-started/advanced-topics/development/)

## Enterprise

Open WebUI Enterprise offers:

- **White-labeling** - Custom branding and theming
- **Dedicated Support** - SLA-backed support
- **Custom Features** - Priority feature development
- **Consulting** - Professional services hours

## Security Policy

### Vulnerability Reporting

Vulnerabilities are reported exclusively through GitHub Security Advisory.

**No private disclosure** - Community-driven security approach.

Report: `https://github.com/open-webui/open-webui/security`

## Community Resources

### Discord Server
Join the community: `https://discord.com/invite/open-webui`

### GitHub Discussions
`https://github.com/open-webui/open-webui/discussions`

### Social Media
- GitHub: `https://github.com/open-webui/open-webui`
- Twitter/X: @openwebui

## Sponsorships

Support Open WebUI development through sponsorships.

**Sponsor Tiers:**
- Individual supporters
- Corporate sponsors
- Enterprise partners

## Roadmap

Future development focuses on:

### Interface Improvements
- Enhanced UX
- Mobile responsiveness
- Custom themes

### Information Retrieval
- Better web search integration
- RAG improvements
- Vector database expansion

### Community
- Better contribution tools
- Translation support
- Accessibility enhancements

---


*This documentation is a comprehensive guide to Open WebUI. For the most up-to-date information, please visit the official documentation at `https://docs.openwebui.com/`*

**Documentation Version:** Complete Scraped Version 2025
**Total Pages Scraped:** 130+
**Last Updated:** 2025-02-07
