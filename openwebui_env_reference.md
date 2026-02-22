# OpenWebUI Environment Variables Comprehensive Reference

> **OpenWebUI Version:** v0.8.3 | **Last Updated:** 2026-02-18

**Related Documents:**
- [Master Reference](OPENWEBUI_MASTER_REFERENCE.md) - Overview and navigation
- [Tools & Functions Guide](openwebui-tools-functions-guide.md) - Development patterns
- [RAG Technical Reference](openwebui_rag_technical_reference.md) - RAG system details
- [Pipelines Guide](OPENWEBUI_PIPELINES_GUIDE.md) - Pipeline deployment
- [API Examples](API_EXAMPLES.md) - Repository-specific API usage

---

> **Version:** v0.8.3 (Latest as of February 2026)
> **Source:** [Official OpenWebUI Documentation](https://docs.openwebui.com/reference/env-configuration/) and [GitHub Source](https://github.com/open-webui/open-webui)

---

## Table of Contents

1. [Important Concepts](#important-concepts)
2. [Core Application Settings](#core-application-settings)
3. [Authentication & Security](#authentication--security)
4. [OAuth / SSO Configuration](#oauth--sso-configuration)
5. [LDAP Configuration](#ldap-configuration)
6. [Database Configuration](#database-configuration)
7. [Redis Configuration](#redis-configuration)
8. [Vector Database (RAG)](#vector-database-rag)
9. [RAG / Knowledge Retrieval](#rag--knowledge-retrieval)
10. [Web Search Configuration](#web-search-configuration)
11. [Model Configuration](#model-configuration)
12. [Image Generation](#image-generation)
13. [Audio / Speech](#audio--speech)
14. [User Permissions](#user-permissions)
15. [Task & Prompt Templates](#task--prompt-templates)
16. [Code Execution](#code-execution)
17. [Storage Providers](#storage-providers)
18. [Logging & Monitoring](#logging--monitoring)
19. [Performance Tuning](#performance-tuning)
20. [WebSocket Configuration](#websocket-configuration)

---

> **See Also:** [RAG Technical Reference](openwebui_rag_technical_reference.md#environment-variables-reference) for RAG-specific environment variables, and [Master Reference](OPENWEBUI_MASTER_REFERENCE.md#3-environment-variables-complete-reference) for overview.

---

## Important Concepts

### PersistentConfig Variables

Environment variables marked as **PersistentConfig** behave specially:

- **First Launch:** All environment variables are used normally
- **After First Launch:** Values are persisted to the database and subsequent restarts use the **database value**, not the environment variable
- **To Override:** Set `ENABLE_PERSISTENT_CONFIG=False` to force using environment variables (not recommended for production)

### Troubleshooting Ignored Environment Variables

If changes to environment variables don't reflect in the UI:

1. **Option 1:** Set `ENABLE_PERSISTENT_CONFIG=False` (temporary, UI changes won't persist)
2. **Option 2:** Update via Admin UI (recommended)
3. **Option 3:** Manual database update: `docker exec -it open-webui sqlite3 /app/backend/data/webui.db "UPDATE config SET data = json_set(data, '$.ENABLE_SIGNUP', json('true'));"`
4. **Option 4:** Reset for fresh install: `docker volume rm open-webui`

---

## Core Application Settings

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENV` | str | `dev` (backend) / `prod` (Docker) | No | Environment mode. Options: `dev` (enables FastAPI docs), `prod` |
| `PORT` | int | `8080` | No | Port to run Open WebUI on |
| `WEBUI_NAME` | str | `Open WebUI` | No | Main WebUI name displayed in the interface |
| `CUSTOM_NAME` | str | `""` | No | Sets WEBUI_NAME but polls api.openwebui.com for metadata |
| `WEBUI_URL` | str | `http://localhost:3000` | **Yes** | Public URL for OAuth/SSO and search engine support. **Must be set before OAuth!** |
| `WEBUI_BUILD_HASH` | str | `dev-build` | No | Git SHA for identifying release version |
| `ENABLE_PERSISTENT_CONFIG` | bool | `True` | No | Controls if database values take precedence over environment variables |
| `DEFAULT_LOCALE` | str | `en` | **Yes** | Default language/locale for the application |
| `ENABLE_EASTER_EGGS` | bool | `True` | No | Enable easter egg features (e.g., "Her" theme) |
| `ENABLE_COMPRESSION_MIDDLEWARE` | bool | `True` | No | Enable gzip compression for HTTP responses |
| `CORS_ALLOW_ORIGIN` | str | `*` | No | CORS allowed origins (semicolon-separated). **Not recommended for production!** |
| `CORS_ALLOW_CUSTOM_SCHEME` | str | `""` | No | Custom URL schemes for CORS (semicolon-separated) |
| `ENABLE_ADMIN_EXPORT` | bool | `True` | No | Allow admins to export data, chats, and database |
| `ENABLE_ADMIN_CHAT_ACCESS` | bool | `True` | No | Allow admins to access user chats |
| `BYPASS_ADMIN_ACCESS_CONTROL` | bool | `True` | No | Admins see all workspace items. Set to `False` to respect access controls |
| `ENABLE_PUBLIC_ACTIVE_USERS_COUNT` | bool | `True` | **Yes** | Show active user count to all users (vs. admins only) |
| `ENABLE_USER_STATUS` | bool | `True` | **Yes** | Enable user status functionality (online/away indicators) |
| `THREAD_POOL_SIZE` | int | `0` (40 threads) | No | FastAPI/AnyIO thread pool size. **Increase for large instances!** |
| `WEBUI_BANNERS` | list | `[]` | **Yes** | JSON array of banner notifications to display |
| `RESPONSE_WATERMARK` | str | `""` | **Yes** | Text added when copying messages (e.g., "AI-generated") |
| `WEBHOOK_URL` | str | `""` | **Yes** | Discord/Slack/Teams webhook URL for notifications |
| `ENABLE_USER_WEBHOOKS` | bool | `True` | **Yes** | Enable user webhooks |
| `ENABLE_CHANNELS` | bool | `False` | **Yes** | Enable channel support |
| `ENABLE_FOLDERS` | bool | `True` | **Yes** | Enable folders feature for organizing chats |
| `FOLDER_MAX_FILE_COUNT` | int | `""` (unlimited) | **Yes** | Max files per folder |
| `ENABLE_NOTES` | bool | `True` | **Yes** | Enable personal notes feature |
| `ENABLE_MEMORIES` | bool | `True` | **Yes** | Enable memory feature for storing user information |
| `ENABLE_SIGNUP` | bool | `True` | **Yes** | Enable user registration |
| `ENABLE_SIGNUP_PASSWORD_CONFIRMATION` | bool | `False` | No | Show password confirmation field on signup |
| `DEFAULT_USER_ROLE` | str | `pending` | **Yes** | Default role for new users: `pending`, `user`, `admin` |
| `DEFAULT_GROUP_ID` | str | `""` | **Yes** | Default group ID for new users |
| `DEFAULT_MODELS` | str | `None` | **Yes** | Default language model ID |
| `DEFAULT_PINNED_MODELS` | str | `""` | **Yes** | Comma-separated list of models to pin for new users |
| `ENABLE_CUSTOM_MODEL_FALLBACK` | bool | `False` | No | Fall back to default model if custom model's base is missing |
| `ENABLE_LOGIN_FORM` | bool | `True` | **Yes** | Show email/password login form |
| `ENABLE_PASSWORD_AUTH` | bool | `True` | No | Enable password authentication. **Only disable when OAuth is configured!** |

### Admin Auto-Creation (Headless Deployment)

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `WEBUI_ADMIN_EMAIL` | str | `""` | Auto-create admin account with this email |
| `WEBUI_ADMIN_PASSWORD` | str | `""` | Password for auto-created admin account |
| `WEBUI_ADMIN_NAME` | str | `Admin` | Display name for auto-created admin |

---

## Authentication & Security

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `WEBUI_AUTH` | bool | `True` | No | Enable authentication (required for security) |
| `WEBUI_SECRET_KEY` | str | `t0p-s3cr3t` | No | JWT signing secret. **Change in production!** |
| `JWT_EXPIRES_IN` | str | `4w` | **Yes** | JWT token expiration. **Security best practice:** Use explicit expiration (e.g., `4w`), not `-1` |
| `ENABLE_API_KEYS` | bool | `False` | **Yes** | Enable API key authentication |
| `ENABLE_API_KEYS_ENDPOINT_RESTRICTIONS` | bool | `False` | **Yes** | Restrict API keys to specific endpoints |
| `API_KEYS_ALLOWED_ENDPOINTS` | str | `""` | **Yes** | Comma-separated list of allowed endpoints |
| `ENABLE_PASSWORD_VALIDATION` | bool | `False` | No | Enable password strength validation |
| `PASSWORD_VALIDATION_REGEX_PATTERN` | regex | `^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)(?=.*[^\w\s]).{8,}$` | No | Regex pattern for password validation |
| `PASSWORD_VALIDATION_HINT` | str | `""` | No | Hint shown for password requirements |
| `BYPASS_MODEL_ACCESS_CONTROL` | bool | `False` | No | All users can access all models regardless of permissions |
| `WEBUI_AUTH_TRUSTED_EMAIL_HEADER` | str | `None` | No | Header for trusted email authentication (reverse proxy) |
| `WEBUI_AUTH_TRUSTED_NAME_HEADER` | str | `None` | No | Header for trusted name authentication |
| `WEBUI_AUTH_TRUSTED_GROUPS_HEADER` | str | `None` | No | Header for trusted groups authentication |
| `WEBUI_AUTH_SIGNOUT_REDIRECT_URL` | str | `None` | No | URL to redirect after signout |
| `WEBUI_SESSION_COOKIE_SAME_SITE` | str | `lax` | No | SameSite attribute for session cookies |
| `WEBUI_SESSION_COOKIE_SECURE` | bool | `False` | No | Secure flag for session cookies |
| `WEBUI_AUTH_COOKIE_SAME_SITE` | str | `lax` | No | SameSite attribute for auth cookies |
| `WEBUI_AUTH_COOKIE_SECURE` | bool | `False` | No | Secure flag for auth cookies |
| `TRUSTED_SIGNATURE_KEY` | str | `""` | No | Key for trusted request signatures |
| `ENABLE_FORWARD_USER_INFO_HEADERS` | bool | `False` | No | Forward user info in HTTP headers to backend |
| `FORWARD_USER_INFO_HEADER_USER_NAME` | str | `X-OpenWebUI-User-Name` | No | Header name for user name |
| `FORWARD_USER_INFO_HEADER_USER_ID` | str | `X-OpenWebUI-User-Id` | No | Header name for user ID |
| `FORWARD_USER_INFO_HEADER_USER_EMAIL` | str | `X-OpenWebUI-User-Email` | No | Header name for user email |
| `FORWARD_USER_INFO_HEADER_USER_ROLE` | str | `X-OpenWebUI-User-Role` | No | Header name for user role |
| `ENABLE_SCIM` | bool | `False` | No | Enable SCIM provisioning |
| `SCIM_TOKEN` | str | `""` | No | SCIM API token |
| `SCIM_AUTH_PROVIDER` | str | `""` | No | OAuth provider name for SCIM externalId |

### SCIM Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_SCIM` | bool | `False` | Enable SCIM user provisioning |
| `SCIM_TOKEN` | str | `""` | Bearer token for SCIM API authentication |
| `SCIM_AUTH_PROVIDER` | str | `""` | OAuth provider name for externalId mapping (e.g., 'microsoft', 'oidc') |

---

## OAuth / SSO Configuration

### General OAuth Settings

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENABLE_OAUTH_SIGNUP` | bool | `False` | **Yes** | Enable OAuth-based signup/login |
| `ENABLE_OAUTH_PERSISTENT_CONFIG` | bool | `False` | No | Allow OAuth settings to persist in database |
| `OAUTH_MERGE_ACCOUNTS_BY_EMAIL` | bool | `False` | **Yes** | Merge existing accounts with same email |
| `ENABLE_OAUTH_ROLE_MANAGEMENT` | bool | `False` | **Yes** | Map OAuth roles to OpenWebUI roles |
| `ENABLE_OAUTH_GROUP_MANAGEMENT` | bool | `False` | **Yes** | Map OAuth groups to OpenWebUI groups |
| `ENABLE_OAUTH_GROUP_CREATION` | bool | `False` | **Yes** | Auto-create groups from OAuth |
| `OAUTH_GROUPS_SEPARATOR` | str | `;` | No | Separator for groups claim |
| `OAUTH_ROLES_SEPARATOR` | str | `,` | No | Separator for roles claim |
| `OAUTH_ALLOWED_DOMAINS` | str | `*` | **Yes** | Comma-separated allowed email domains |
| `OAUTH_ALLOWED_ROLES` | str | `user,admin` | **Yes** | Comma-separated allowed roles |
| `OAUTH_ADMIN_ROLES` | str | `admin` | **Yes** | Comma-separated roles that get admin access |
| `OAUTH_ROLES_CLAIM` | str | `roles` | **Yes** | JWT claim name for roles |
| `OAUTH_GROUPS_CLAIM` | str | `groups` | **Yes** | JWT claim name for groups |
| `OAUTH_BLOCKED_GROUPS` | str | `[]` | **Yes** | JSON array of blocked group names |
| `OAUTH_UPDATE_PICTURE_ON_LOGIN` | bool | `False` | **Yes** | Update user picture from OAuth on each login |
| `OAUTH_ACCESS_TOKEN_REQUEST_INCLUDE_CLIENT_ID` | bool | `False` | No | Include client_id in token request |
| `OAUTH_AUDIENCE` | str | `""` | **Yes** | OAuth audience parameter |
| `ENABLE_OAUTH_EMAIL_FALLBACK` | bool | `False` | No | Allow email fallback for OAuth |
| `ENABLE_OAUTH_ID_TOKEN_COOKIE` | bool | `True` | No | Store ID token in cookie |
| `OAUTH_CLIENT_INFO_ENCRYPTION_KEY` | str | `WEBUI_SECRET_KEY` | No | Key for encrypting OAuth client info |
| `OAUTH_SESSION_TOKEN_ENCRYPTION_KEY` | str | `WEBUI_SECRET_KEY` | No | Key for encrypting session tokens |
| `ENABLE_OAUTH_TOKEN_EXCHANGE` | bool | `False` | No | Enable token exchange for external apps |

### Google OAuth

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `GOOGLE_CLIENT_ID` | str | `""` | **Yes** |
| `GOOGLE_CLIENT_SECRET` | str | `""` | **Yes** |
| `GOOGLE_OAUTH_SCOPE` | str | `openid email profile` | **Yes** |
| `GOOGLE_REDIRECT_URI` | str | `""` | **Yes** |

### Microsoft OAuth

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `MICROSOFT_CLIENT_ID` | str | `""` | **Yes** |
| `MICROSOFT_CLIENT_SECRET` | str | `""` | **Yes** |
| `MICROSOFT_CLIENT_TENANT_ID` | str | `""` | **Yes** |
| `MICROSOFT_CLIENT_LOGIN_BASE_URL` | str | `https://login.microsoftonline.com` | **Yes** |
| `MICROSOFT_CLIENT_PICTURE_URL` | str | `https://graph.microsoft.com/v1.0/me/photo/$value` | **Yes** |
| `MICROSOFT_OAUTH_SCOPE` | str | `openid email profile` | **Yes** |
| `MICROSOFT_REDIRECT_URI` | str | `""` | **Yes** |

### GitHub OAuth

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `GITHUB_CLIENT_ID` | str | `""` | **Yes** |
| `GITHUB_CLIENT_SECRET` | str | `""` | **Yes** |
| `GITHUB_CLIENT_SCOPE` | str | `user:email` | **Yes** |
| `GITHUB_CLIENT_REDIRECT_URI` | str | `""` | **Yes** |

### Generic OIDC OAuth

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `OAUTH_CLIENT_ID` | str | `""` | **Yes** | OIDC client ID |
| `OAUTH_CLIENT_SECRET` | str | `""` | **Yes** | OIDC client secret |
| `OPENID_PROVIDER_URL` | str | `""` | **Yes** | OIDC discovery URL |
| `OPENID_REDIRECT_URI` | str | `""` | **Yes** | Redirect URI |
| `OAUTH_SCOPES` | str | `openid email profile` | **Yes** | Requested scopes |
| `OAUTH_TIMEOUT` | int | `""` | **Yes** | OAuth request timeout |
| `OAUTH_TOKEN_ENDPOINT_AUTH_METHOD` | str | `None` | **Yes** | Token endpoint auth method |
| `OAUTH_CODE_CHALLENGE_METHOD` | str | `None` | **Yes** | PKCE code challenge method (S256) |
| `OAUTH_PROVIDER_NAME` | str | `SSO` | **Yes** | Display name for provider |
| `OAUTH_SUB_CLAIM` | str | `None` | **Yes** | Claim for subject ID |
| `OAUTH_USERNAME_CLAIM` | str | `name` | **Yes** | Claim for username |
| `OAUTH_PICTURE_CLAIM` | str | `picture` | **Yes** | Claim for avatar URL |
| `OAUTH_EMAIL_CLAIM` | str | `email` | **Yes** | Claim for email |

### Feishu OAuth

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `FEISHU_CLIENT_ID` | str | `""` | **Yes** |
| `FEISHU_CLIENT_SECRET` | str | `""` | **Yes** |
| `FEISHU_OAUTH_SCOPE` | str | `contact:user.base:readonly` | **Yes** |
| `FEISHU_REDIRECT_URI` | str | `""` | **Yes** |

---

## LDAP Configuration

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENABLE_LDAP` | bool | `False` | **Yes** | Enable LDAP authentication |
| `LDAP_SERVER_LABEL` | str | `LDAP Server` | **Yes** | Display label for LDAP server |
| `LDAP_SERVER_HOST` | str | `localhost` | **Yes** | LDAP server hostname |
| `LDAP_SERVER_PORT` | int | `389` | **Yes** | LDAP server port |
| `LDAP_ATTRIBUTE_FOR_MAIL` | str | `mail` | **Yes** | LDAP attribute for email |
| `LDAP_ATTRIBUTE_FOR_USERNAME` | str | `uid` | **Yes** | LDAP attribute for username |
| `LDAP_APP_DN` | str | `""` | **Yes** | DN for LDAP bind user |
| `LDAP_APP_PASSWORD` | str | `""` | **Yes** | Password for LDAP bind user |
| `LDAP_SEARCH_BASE` | str | `""` | **Yes** | Base DN for user search |
| `LDAP_SEARCH_FILTERS` | str | `""` | **Yes** | Additional LDAP search filters |
| `LDAP_USE_TLS` | bool | `True` | **Yes** | Use TLS for LDAP connections |
| `LDAP_CA_CERT_FILE` | str | `""` | **Yes** | Path to CA certificate |
| `LDAP_VALIDATE_CERT` | bool | `True` | **Yes** | Validate LDAP server certificate |
| `LDAP_CIPHERS` | str | `ALL` | **Yes** | Allowed SSL ciphers |
| `ENABLE_LDAP_GROUP_MANAGEMENT` | bool | `False` | **Yes** | Enable LDAP group mapping |
| `ENABLE_LDAP_GROUP_CREATION` | bool | `False` | **Yes** | Auto-create groups from LDAP |
| `LDAP_ATTRIBUTE_FOR_GROUPS` | str | `memberOf` | **Yes** | LDAP attribute for group membership |

---

## Database Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `DATABASE_URL` | str | `sqlite:///{DATA_DIR}/webui.db` | Full database connection URL |
| `DATABASE_TYPE` | str | `None` | Database type (postgres, mysql, etc.) |
| `DATABASE_USER` | str | `None` | Database username |
| `DATABASE_PASSWORD` | str | `None` | Database password |
| `DATABASE_HOST` | str | `None` | Database host |
| `DATABASE_PORT` | str | `None` | Database port |
| `DATABASE_NAME` | str | `None` | Database name |
| `DATABASE_SCHEMA` | str | `None` | Database schema (PostgreSQL) |
| `DATABASE_POOL_SIZE` | int | `None` | Connection pool size |
| `DATABASE_POOL_MAX_OVERFLOW` | int | `0` | Max overflow connections |
| `DATABASE_POOL_TIMEOUT` | int | `30` | Pool connection timeout (seconds) |
| `DATABASE_POOL_RECYCLE` | int | `3600` | Connection recycle time (seconds) |
| `DATABASE_ENABLE_SQLITE_WAL` | bool | `False` | Enable SQLite WAL mode |
| `DATABASE_ENABLE_SESSION_SHARING` | bool | `False` | Reuse existing DB sessions |
| `DATABASE_USER_ACTIVE_STATUS_UPDATE_INTERVAL` | float | `None` | Interval for updating user active status |
| `ENABLE_DB_MIGRATIONS` | bool | `True` | Run database migrations on startup |

---

## Redis Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `REDIS_URL` | str | `""` | Redis connection URL (required for multi-instance) |
| `REDIS_CLUSTER` | bool | `False` | Use Redis Cluster mode |
| `REDIS_KEY_PREFIX` | str | `open-webui` | Prefix for Redis keys |
| `REDIS_SENTINEL_HOSTS` | str | `""` | Comma-separated Sentinel hosts |
| `REDIS_SENTINEL_PORT` | str | `26379` | Sentinel port |
| `REDIS_SENTINEL_MAX_RETRY_COUNT` | int | `2` | Max retries for Sentinel failover |
| `REDIS_SOCKET_CONNECT_TIMEOUT` | float | `""` | Socket connect timeout |
| `REDIS_RECONNECT_DELAY` | float | `None` | Delay between reconnection attempts |

---

## Vector Database (RAG)

### Vector DB Selection

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `VECTOR_DB` | str | `chroma` | Vector DB: `chroma`, `milvus`, `qdrant`, `weaviate`, `opensearch`, `elasticsearch`, `pgvector`, `pinecone`, `oracle23ai` |

### Chroma Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `CHROMA_DATA_PATH` | str | `{DATA_DIR}/vector_db` | Chroma data directory |
| `CHROMA_TENANT` | str | `default_tenant` | Chroma tenant |
| `CHROMA_DATABASE` | str | `default_database` | Chroma database |
| `CHROMA_HTTP_HOST` | str | `""` | Chroma HTTP host (for client mode) |
| `CHROMA_HTTP_PORT` | int | `8000` | Chroma HTTP port |
| `CHROMA_CLIENT_AUTH_PROVIDER` | str | `""` | Auth provider |
| `CHROMA_CLIENT_AUTH_CREDENTIALS` | str | `""` | Auth credentials |
| `CHROMA_HTTP_HEADERS` | str | `""` | HTTP headers (format: `header1=value1,header2=value2`) |
| `CHROMA_HTTP_SSL` | bool | `False` | Use SSL for HTTP |

### Milvus Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `MILVUS_URI` | str | `{DATA_DIR}/vector_db/milvus.db` | Milvus connection URI |
| `MILVUS_DB` | str | `default` | Milvus database name |
| `MILVUS_TOKEN` | str | `None` | Milvus authentication token |
| `MILVUS_INDEX_TYPE` | str | `HNSW` | Index type: `HNSW`, `IVF_FLAT`, `DISKANN` |
| `MILVUS_METRIC_TYPE` | str | `COSINE` | Metric type: `COSINE`, `L2`, `IP` |
| `MILVUS_HNSW_M` | int | `16` | HNSW M parameter |
| `MILVUS_HNSW_EFCONSTRUCTION` | int | `100` | HNSW efConstruction parameter |
| `MILVUS_IVF_FLAT_NLIST` | int | `128` | IVF Flat nlist parameter |
| `MILVUS_DISKANN_MAX_DEGREE` | int | `56` | DiskANN max degree |
| `MILVUS_DISKANN_SEARCH_LIST_SIZE` | int | `100` | DiskANN search list size |
| `ENABLE_MILVUS_MULTITENANCY_MODE` | bool | `False` | Enable multi-tenancy |
| `MILVUS_COLLECTION_PREFIX` | str | `open_webui` | Collection name prefix |

### Qdrant Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `QDRANT_URI` | str | `None` | Qdrant server URI |
| `QDRANT_API_KEY` | str | `None` | Qdrant API key |
| `QDRANT_ON_DISK` | bool | `False` | Store vectors on disk |
| `QDRANT_PREFER_GRPC` | bool | `False` | Prefer gRPC over HTTP |
| `QDRANT_GRPC_PORT` | int | `6334` | gRPC port |
| `QDRANT_TIMEOUT` | int | `5` | Request timeout (seconds) |
| `QDRANT_HNSW_M` | int | `16` | HNSW M parameter |
| `ENABLE_QDRANT_MULTITENANCY_MODE` | bool | `True` | Enable multi-tenancy |
| `QDRANT_COLLECTION_PREFIX` | str | `open-webui` | Collection name prefix |

### Weaviate Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `WEAVIATE_HTTP_HOST` | str | `""` | HTTP host |
| `WEAVIATE_GRPC_HOST` | str | `""` | gRPC host |
| `WEAVIATE_HTTP_PORT` | int | `8080` | HTTP port |
| `WEAVIATE_GRPC_PORT` | int | `50051` | gRPC port |
| `WEAVIATE_API_KEY` | str | `None` | API key |
| `WEAVIATE_HTTP_SECURE` | bool | `False` | Use HTTPS |
| `WEAVIATE_GRPC_SECURE` | bool | `False` | Use gRPC TLS |
| `WEAVIATE_SKIP_INIT_CHECKS` | bool | `False` | Skip initialization checks |

### OpenSearch Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `OPENSEARCH_URI` | str | `https://localhost:9200` | OpenSearch URI |
| `OPENSEARCH_SSL` | bool | `True` | Use SSL |
| `OPENSEARCH_CERT_VERIFY` | bool | `False` | Verify certificates |
| `OPENSEARCH_USERNAME` | str | `None` | Username |
| `OPENSEARCH_PASSWORD` | str | `None` | Password |

### Elasticsearch Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ELASTICSEARCH_URL` | str | `https://localhost:9200` | Elasticsearch URL |
| `ELASTICSEARCH_CA_CERTS` | str | `None` | Path to CA certificates |
| `ELASTICSEARCH_API_KEY` | str | `None` | API key |
| `ELASTICSEARCH_USERNAME` | str | `None` | Username |
| `ELASTICSEARCH_PASSWORD` | str | `None` | Password |
| `ELASTICSEARCH_CLOUD_ID` | str | `None` | Cloud ID |
| `SSL_ASSERT_FINGERPRINT` | str | `None` | SSL fingerprint |
| `ELASTICSEARCH_INDEX_PREFIX` | str | `open_webui_collections` | Index prefix |

### Pgvector Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `PGVECTOR_DB_URL` | str | `DATABASE_URL` | PostgreSQL connection URL |
| `PGVECTOR_INITIALIZE_MAX_VECTOR_LENGTH` | int | `1536` | Max vector dimension |
| `PGVECTOR_USE_HALFVEC` | bool | `False` | Use halfvec for >2000 dims |
| `PGVECTOR_CREATE_EXTENSION` | bool | `True` | Auto-create pgvector extension |
| `PGVECTOR_PGCRYPTO` | bool | `False` | Enable pgcrypto encryption |
| `PGVECTOR_PGCRYPTO_KEY` | str | `None` | pgcrypto encryption key |
| `PGVECTOR_POOL_SIZE` | int | `None` | Connection pool size |
| `PGVECTOR_POOL_MAX_OVERFLOW` | int | `0` | Max overflow connections |
| `PGVECTOR_POOL_TIMEOUT` | int | `30` | Pool timeout |
| `PGVECTOR_POOL_RECYCLE` | int | `3600` | Connection recycle time |
| `PGVECTOR_INDEX_METHOD` | str | `""` | Index method: `ivfflat`, `hnsw` |
| `PGVECTOR_HNSW_M` | int | `16` | HNSW M parameter |
| `PGVECTOR_HNSW_EF_CONSTRUCTION` | int | `64` | HNSW ef_construction |
| `PGVECTOR_IVFFLAT_LISTS` | int | `100` | IVFFlat lists |

### Pinecone Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `PINECONE_API_KEY` | str | `None` | Pinecone API key |
| `PINECONE_ENVIRONMENT` | str | `None` | Pinecone environment |
| `PINECONE_INDEX_NAME` | str | `open-webui-index` | Index name |
| `PINECONE_DIMENSION` | int | `1536` | Vector dimension |
| `PINECONE_METRIC` | str | `cosine` | Distance metric |
| `PINECONE_CLOUD` | str | `aws` | Cloud provider: `aws`, `gcp`, `azure` |

### Oracle23ai Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ORACLE_DB_USE_WALLET` | bool | `False` | Use wallet authentication |
| `ORACLE_DB_USER` | str | `None` | Database user |
| `ORACLE_DB_PASSWORD` | str | `None` | Database password |
| `ORACLE_DB_DSN` | str | `None` | Database DSN |
| `ORACLE_WALLET_DIR` | str | `None` | Wallet directory |
| `ORACLE_WALLET_PASSWORD` | str | `None` | Wallet password |
| `ORACLE_VECTOR_LENGTH` | int | `768` | Vector dimension |
| `ORACLE_DB_POOL_MIN` | int | `2` | Min pool size |
| `ORACLE_DB_POOL_MAX` | int | `10` | Max pool size |
| `ORACLE_DB_POOL_INCREMENT` | int | `1` | Pool increment |

### S3 Vector Store

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `S3_VECTOR_BUCKET_NAME` | str | `None` | S3 bucket for vectors |
| `S3_VECTOR_REGION` | str | `None` | S3 region |

---

## RAG / Knowledge Retrieval

### Core RAG Settings

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `BYPASS_EMBEDDING_AND_RETRIEVAL` | bool | `False` | **Yes** | Skip embedding/retrieval (use full document) |
| `RAG_TOP_K` | int | `5` | **Yes** | Number of chunks to retrieve |
| `RAG_TOP_K_RERANKER` | int | `3` | **Yes** | Number of chunks after reranking |
| `RAG_RELEVANCE_THRESHOLD` | float | `0.0` | **Yes** | Minimum relevance score threshold |
| `RAG_HYBRID_BM25_WEIGHT` | float | `0.5` | **Yes** | Weight for BM25 in hybrid search |
| `ENABLE_RAG_HYBRID_SEARCH` | bool | `False` | **Yes** | Enable hybrid (dense + sparse) search |
| `ENABLE_RAG_HYBRID_SEARCH_ENRICHED_TEXTS` | bool | `False` | **Yes** | Use enriched texts for hybrid search |
| `RAG_FULL_CONTEXT` | bool | `False` | **Yes** | Use full document context (no chunking) |
| `RAG_FILE_MAX_COUNT` | int | `None` | **Yes** | Max files per upload |
| `RAG_FILE_MAX_SIZE` | int | `None` | **Yes** | Max file size in bytes |
| `RAG_ALLOWED_FILE_EXTENSIONS` | str | `""` | **Yes** | Comma-separated allowed extensions |
| `PDF_EXTRACT_IMAGES` | bool | `False` | **Yes** | Extract images from PDFs |
| `PDF_LOADER_MODE` | str | `page` | **Yes** | PDF loading mode: `page`, `single` |
| `RAG_EMBEDDING_ENGINE` | str | `""` | **Yes** | Embedding engine: `ollama`, `openai`, `hf` |
| `RAG_EMBEDDING_MODEL` | str | `sentence-transformers/all-MiniLM-L6-v2` | **Yes** | Embedding model name |
| `RAG_SYSTEM_CONTEXT` | bool | `true` | **Yes** | Inject RAG context into system prompt for KV cache optimization |
| `RAG_EMBEDDING_MODEL_AUTO_UPDATE` | bool | `True` (if not OFFLINE) | No | Auto-update embedding model |
| `RAG_EMBEDDING_MODEL_TRUST_REMOTE_CODE` | bool | `True` | No | Trust remote code for embeddings |
| `RAG_EMBEDDING_BATCH_SIZE` | int | `1` | **Yes** | Batch size for embeddings |
| `ENABLE_ASYNC_EMBEDDING` | bool | `True` | **Yes** | Enable async embedding generation |
| `RAG_EMBEDDING_QUERY_PREFIX` | str | `None` | No | Prefix for query embeddings |
| `RAG_EMBEDDING_CONTENT_PREFIX` | str | `None` | No | Prefix for content embeddings |
| `RAG_EMBEDDING_PREFIX_FIELD_NAME` | str | `None` | No | Field name for prefix |
| `RAG_EMBEDDING_TIMEOUT` | int | `None` | No | Embedding request timeout |

### Reranking

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `RAG_RERANKING_ENGINE` | str | `""` | **Yes** | Reranking engine: `hf`, `openai` |
| `RAG_RERANKING_MODEL` | str | `""` | **Yes** | Reranking model name |
| `RAG_RERANKING_MODEL_AUTO_UPDATE` | bool | `True` (if not OFFLINE) | No | Auto-update reranking model |
| `RAG_RERANKING_MODEL_TRUST_REMOTE_CODE` | bool | `True` | No | Trust remote code |
| `RAG_EXTERNAL_RERANKER_URL` | str | `""` | **Yes** | External reranker API URL |
| `RAG_EXTERNAL_RERANKER_API_KEY` | str | `""` | **Yes** | External reranker API key |
| `RAG_EXTERNAL_RERANKER_TIMEOUT` | str | `""` | **Yes** | External reranker timeout |

### Text Splitting

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `RAG_TEXT_SPLITTER` | str | `""` | **Yes** | Text splitter: `token`, `character`, `tiktoken` |
| `CHUNK_SIZE` | int | `1500` | **Yes** | Chunk size in characters/tokens |
| `CHUNK_MIN_SIZE_TARGET` | int | `0` | **Yes** | Minimum chunk size target |
| `CHUNK_OVERLAP` | int | `150` | **Yes** | Overlap between chunks |
| `ENABLE_MARKDOWN_HEADER_TEXT_SPLITTER` | bool | `True` | **Yes** | Split by markdown headers |
| `TIKTOKEN_ENCODING_NAME` | str | `cl100k_base` | **Yes** | Tiktoken encoding |
| `TIKTOKEN_CACHE_DIR` | str | `{CACHE_DIR}/tiktoken` | No | Tiktoken cache directory |

### RAG Template

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `RAG_TEMPLATE` | str | *(see docs)* | **Yes** | System prompt template for RAG |

### External Document Loaders

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `CONTENT_EXTRACTION_ENGINE` | str | `""` | **Yes** | Engine: `tika`, `docling`, `marker`, `mineru`, `mistral` |
| `TIKA_SERVER_URL` | str | `http://tika:9998/tika` | **Yes** | Tika server URL |
| `DOCLING_SERVER_URL` | str | `http://docling:5001` | **Yes** | Docling server URL |
| `DOCLING_API_KEY` | str | `""` | **Yes** | Docling API key |
| `DOCLING_PARAMS` | str | `{}` | **Yes** | JSON parameters for Docling |
| `DOCUMENT_INTELLIGENCE_ENDPOINT` | str | `""` | **Yes** | Azure Document Intelligence endpoint |
| `DOCUMENT_INTELLIGENCE_KEY` | str | `""` | **Yes** | Azure Document Intelligence key |
| `DOCUMENT_INTELLIGENCE_MODEL` | str | `prebuilt-layout` | **Yes** | Model to use |
| `MISTRAL_OCR_API_BASE_URL` | str | `https://api.mistral.ai/v1` | **Yes** | Mistral OCR API URL |
| `MISTRAL_OCR_API_KEY` | str | `""` | **Yes** | Mistral OCR API key |

### Marker Configuration

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `DATALAB_MARKER_API_KEY` | str | `""` | **Yes** |
| `DATALAB_MARKER_API_BASE_URL` | str | `""` | **Yes** |
| `DATALAB_MARKER_ADDITIONAL_CONFIG` | str | `""` | **Yes** |
| `DATALAB_MARKER_USE_LLM` | bool | `False` | **Yes** |
| `DATALAB_MARKER_SKIP_CACHE` | bool | `False` | **Yes** |
| `DATALAB_MARKER_FORCE_OCR` | bool | `False` | **Yes** |
| `DATALAB_MARKER_PAGINATE` | bool | `False` | **Yes** |
| `DATALAB_MARKER_STRIP_EXISTING_OCR` | bool | `False` | **Yes** |
| `DATALAB_MARKER_DISABLE_IMAGE_EXTRACTION` | bool | `False` | **Yes** |
| `DATALAB_MARKER_FORMAT_LINES` | bool | `False` | **Yes** |
| `DATALAB_MARKER_OUTPUT_FORMAT` | str | `markdown` | **Yes** |

### MinerU Configuration

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `MINERU_API_MODE` | str | `local` | **Yes** |
| `MINERU_API_URL` | str | `http://localhost:8000` | **Yes** |
| `MINERU_API_TIMEOUT` | str | `300` | **Yes** |
| `MINERU_API_KEY` | str | `""` | **Yes** |
| `MINERU_PARAMS` | str | `{}` | **Yes** |

### External Loaders

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `EXTERNAL_DOCUMENT_LOADER_URL` | str | `""` | **Yes** |
| `EXTERNAL_DOCUMENT_LOADER_API_KEY` | str | `""` | **Yes** |

### RAG API Settings

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `RAG_OPENAI_API_BASE_URL` | str | `OPENAI_API_BASE_URL` | **Yes** |
| `RAG_OPENAI_API_KEY` | str | `OPENAI_API_KEY` | **Yes** |
| `RAG_AZURE_OPENAI_BASE_URL` | str | `""` | **Yes** |
| `RAG_AZURE_OPENAI_API_KEY` | str | `""` | **Yes** |
| `RAG_AZURE_OPENAI_API_VERSION` | str | `""` | **Yes** |
| `RAG_OLLAMA_BASE_URL` | str | `OLLAMA_BASE_URL` | **Yes** |
| `RAG_OLLAMA_API_KEY` | str | `""` | **Yes** |

### YouTube Loader

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `YOUTUBE_LOADER_LANGUAGE` | str | `en` | **Yes** |
| `YOUTUBE_LOADER_PROXY_URL` | str | `""` | **Yes** |

### Web Fetch

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_RAG_LOCAL_WEB_FETCH` | bool | `False` | Allow fetching local URLs |
| `WEB_FETCH_FILTER_LIST` | str | `*cloud metadata*` | Comma-separated URL filter patterns |

---

## Web Search Configuration

### General Web Search

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENABLE_WEB_SEARCH` | bool | `False` | **Yes** | Enable web search feature |
| `WEB_SEARCH_ENGINE` | str | `""` | **Yes** | Engine: `searxng`, `google_pse`, `brave`, `kagi`, `mojeek`, `bocha`, `serpstack`, `serper`, `serply`, `searchapi`, `serpapi`, `bing`, `azure_ai_search`, `exa`, `perplexity`, `tavily`, `jina`, `ddgs`, `yandex` |
| `BYPASS_WEB_SEARCH_EMBEDDING_AND_RETRIEVAL` | bool | `False` | **Yes** | Skip embedding for web search results |
| `BYPASS_WEB_SEARCH_WEB_LOADER` | bool | `False` | **Yes** | Skip web page loading |
| `WEB_SEARCH_RESULT_COUNT` | int | `3` | **Yes** | Number of search results |
| `WEB_SEARCH_DOMAIN_FILTER_LIST` | list | `[]` | **Yes** | Filter results to these domains |
| `WEB_SEARCH_CONCURRENT_REQUESTS` | int | `0` | **Yes** | Concurrent web requests (0=unlimited) |
| `WEB_SEARCH_TRUST_ENV` | bool | `False` | **Yes** | Trust environment for proxy settings |

### Web Loader

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `WEB_LOADER_ENGINE` | str | `""` | **Yes** | Loader: `playwright`, `firecrawl`, `external` |
| `WEB_LOADER_CONCURRENT_REQUESTS` | int | `10` | **Yes** | Concurrent loader requests |
| `WEB_LOADER_TIMEOUT` | str | `""` | **Yes** | Loader timeout |
| `ENABLE_WEB_LOADER_SSL_VERIFICATION` | bool | `True` | **Yes** | Verify SSL certificates |

### Search Engine APIs

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `SEARXNG_QUERY_URL` | str | `""` | **Yes** |
| `SEARXNG_LANGUAGE` | str | `all` | **Yes** |
| `YACY_QUERY_URL` | str | `""` | **Yes** |
| `YACY_USERNAME` | str | `""` | **Yes** |
| `YACY_PASSWORD` | str | `""` | **Yes** |
| `GOOGLE_PSE_API_KEY` | str | `""` | **Yes** |
| `GOOGLE_PSE_ENGINE_ID` | str | `""` | **Yes** |
| `BRAVE_SEARCH_API_KEY` | str | `""` | **Yes** |
| `KAGI_SEARCH_API_KEY` | str | `""` | **Yes** |
| `MOJEEK_SEARCH_API_KEY` | str | `""` | **Yes** |
| `BOCHA_SEARCH_API_KEY` | str | `""` | **Yes** |
| `SERPSTACK_API_KEY` | str | `""` | **Yes** |
| `SERPSTACK_HTTPS` | bool | `True` | **Yes** |
| `SERPER_API_KEY` | str | `""` | **Yes** |
| `SERPLY_API_KEY` | str | `""` | **Yes** |
| `DDGS_BACKEND` | str | `auto` | **Yes** |
| `JINA_API_KEY` | str | `""` | **Yes** |
| `JINA_API_BASE_URL` | str | `""` | **Yes** |
| `SEARCHAPI_API_KEY` | str | `""` | **Yes** |
| `SEARCHAPI_ENGINE` | str | `""` | **Yes** |
| `SERPAPI_API_KEY` | str | `""` | **Yes** |
| `SERPAPI_ENGINE` | str | `""` | **Yes** |
| `BING_SEARCH_V7_ENDPOINT` | str | `https://api.bing.microsoft.com/v7.0/search` | **Yes** |
| `BING_SEARCH_V7_SUBSCRIPTION_KEY` | str | `""` | **Yes** |
| `AZURE_AI_SEARCH_API_KEY` | str | `""` | **Yes** |
| `AZURE_AI_SEARCH_ENDPOINT` | str | `""` | **Yes** |
| `AZURE_AI_SEARCH_INDEX_NAME` | str | `""` | **Yes** |
| `EXA_API_KEY` | str | `""` | **Yes** |
| `PERPLEXITY_API_KEY` | str | `""` | **Yes** |
| `PERPLEXITY_MODEL` | str | `sonar` | **Yes** |
| `PERPLEXITY_SEARCH_CONTEXT_USAGE` | str | `medium` | **Yes** |
| `PERPLEXITY_SEARCH_API_URL` | str | `https://api.perplexity.ai/search` | **Yes** |
| `SOUGOU_API_SID` | str | `""` | **Yes** |
| `SOUGOU_API_SK` | str | `""` | **Yes** |
| `TAVILY_API_KEY` | str | `""` | **Yes** |
| `TAVILY_EXTRACT_DEPTH` | str | `basic` | **Yes** |

### External Web Search/Loader

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `EXTERNAL_WEB_SEARCH_URL` | str | `""` | **Yes** |
| `EXTERNAL_WEB_SEARCH_API_KEY` | str | `""` | **Yes** |
| `EXTERNAL_WEB_LOADER_URL` | str | `""` | **Yes** |
| `EXTERNAL_WEB_LOADER_API_KEY` | str | `""` | **Yes** |

### Playwright

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `PLAYWRIGHT_WS_URL` | str | `""` | **Yes** |
| `PLAYWRIGHT_TIMEOUT` | int | `10000` | **Yes** |

### Firecrawl

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `FIRECRAWL_API_KEY` | str | `""` | **Yes** |
| `FIRECRAWL_API_BASE_URL` | str | `https://api.firecrawl.dev` | **Yes** |
| `FIRECRAWL_TIMEOUT` | str | `""` | **Yes** |

### Yandex

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `YANDEX_WEB_SEARCH_URL` | str | `""` | **Yes** |
| `YANDEX_WEB_SEARCH_API_KEY` | str | `""` | **Yes** |
| `YANDEX_WEB_SEARCH_CONFIG` | str | `""` | **Yes** |

---

## Model Configuration

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENABLE_OLLAMA_API` | bool | `True` | **Yes** | Enable Ollama API |
| `OLLAMA_BASE_URL` | str | `http://localhost:11434` | No | Ollama server URL |
| `OLLAMA_BASE_URLS` | str | `OLLAMA_BASE_URL` | **Yes** | Semicolon-separated Ollama URLs for load balancing |
| `OLLAMA_API_CONFIGS` | dict | `{}` | **Yes** | Additional Ollama API configurations |
| `K8S_FLAG` | str | `""` | No | Kubernetes deployment flag |
| `USE_OLLAMA_DOCKER` | str | `false` | No | Use bundled Ollama in Docker |
| `ENABLE_OPENAI_API` | bool | `True` | **Yes** | Enable OpenAI-compatible API |
| `OPENAI_API_KEY` | str | `""` | **Yes** | OpenAI API key |
| `OPENAI_API_KEYS` | str | `OPENAI_API_KEY` | **Yes** | Semicolon-separated keys for load balancing |
| `OPENAI_API_BASE_URL` | str | `https://api.openai.com/v1` | **Yes** | OpenAI API base URL |
| `OPENAI_API_BASE_URLS` | str | `OPENAI_API_BASE_URL` | **Yes** | Semicolon-separated URLs |
| `OPENAI_API_CONFIGS` | dict | `{}` | **Yes** | Additional API configurations |
| `GEMINI_API_KEY` | str | `""` | No | Google Gemini API key |
| `GEMINI_API_BASE_URL` | str | `""` | No | Gemini API base URL |
| `ENABLE_BASE_MODELS_CACHE` | bool | `False` | **Yes** | Cache base models list |
| `MODELS_CACHE_TTL` | int | `1` | No | Cache TTL in seconds |
| `ENABLE_DIRECT_CONNECTIONS` | bool | `False` | **Yes** | Enable direct model connections |
| `MODEL_ORDER_LIST` | list | `[]` | **Yes** | Custom model ordering |

---

## Image Generation

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENABLE_IMAGE_GENERATION` | bool | `False` | **Yes** | Enable image generation |
| `IMAGE_GENERATION_ENGINE` | str | `openai` | **Yes** | Engine: `openai`, `automatic1111`, `comfyui`, `gemini` |
| `IMAGE_GENERATION_MODEL` | str | `""` | **Yes** | Model name |
| `IMAGE_SIZE` | str | `512x512` | **Yes** | Generated image size |
| `IMAGE_STEPS` | int | `50` | **Yes** | Generation steps |
| `ENABLE_IMAGE_PROMPT_GENERATION` | bool | `True` | **Yes** | Auto-generate image prompts |
| `IMAGE_AUTO_SIZE_MODELS_REGEX_PATTERN` | str | `^gpt-image` | No | Regex for auto-size models |
| `IMAGE_URL_RESPONSE_MODELS_REGEX_PATTERN` | str | `^gpt-image` | No | Regex for URL response models |
| `ENABLE_IMAGE_EDIT` | bool | `False` | **Yes** | Enable image editing |
| `IMAGE_EDIT_ENGINE` | str | `openai` | **Yes** | Edit engine |
| `IMAGE_EDIT_MODEL` | str | `""` | **Yes** | Edit model |
| `IMAGE_EDIT_SIZE` | str | `""` | **Yes** | Edit output size |

### Automatic1111

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `AUTOMATIC1111_BASE_URL` | str | `""` | **Yes** |
| `AUTOMATIC1111_API_AUTH` | str | `""` | **Yes** |
| `AUTOMATIC1111_PARAMS` | str | `{}` | **Yes** |

### ComfyUI

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `COMFYUI_BASE_URL` | str | `""` | **Yes** |
| `COMFYUI_API_KEY` | str | `""` | **Yes** |
| `COMFYUI_WORKFLOW` | str | *(default)* | **Yes** |
| `COMFYUI_WORKFLOW_NODES` | str | `[]` | **Yes** |

### OpenAI Images

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `IMAGES_OPENAI_API_BASE_URL` | str | `OPENAI_API_BASE_URL` | **Yes** |
| `IMAGES_OPENAI_API_VERSION` | str | `""` | **Yes** |
| `IMAGES_OPENAI_API_KEY` | str | `OPENAI_API_KEY` | **Yes** |
| `IMAGES_OPENAI_API_PARAMS` | str | `{}` | **Yes** |

### Gemini Images

| Variable | Type | Default | Persistent |
|----------|------|---------|------------|
| `IMAGES_GEMINI_API_BASE_URL` | str | `GEMINI_API_BASE_URL` | **Yes** |
| `IMAGES_GEMINI_API_KEY` | str | `GEMINI_API_KEY` | **Yes** |
| `IMAGES_GEMINI_ENDPOINT_METHOD` | str | `""` | **Yes** |

---

## Audio / Speech

### Speech-to-Text (STT)

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `AUDIO_STT_ENGINE` | str | `""` | **Yes** | STT engine: `whisper`, `openai`, `deepgram`, `azure`, `mistral` |
| `AUDIO_STT_MODEL` | str | `""` | **Yes** | STT model name |
| `WHISPER_MODEL` | str | `base` | **Yes** | Whisper model size |
| `WHISPER_COMPUTE_TYPE` | str | `int8` | No | Compute type: `int8`, `float16`, `float32` |
| `WHISPER_MODEL_DIR` | str | `{CACHE_DIR}/whisper/models` | No | Whisper model directory |
| `WHISPER_MODEL_AUTO_UPDATE` | bool | `False` | No | Auto-update whisper model |
| `WHISPER_VAD_FILTER` | bool | `False` | No | Enable voice activity detection |
| `WHISPER_MULTILINGUAL` | bool | `False` | No | Enable multilingual support |
| `WHISPER_LANGUAGE` | str | `""` | No | Force specific language |
| `DEEPGRAM_API_KEY` | str | `""` | **Yes** | Deepgram API key |
| `AUDIO_STT_OPENAI_API_BASE_URL` | str | `OPENAI_API_BASE_URL` | **Yes** |
| `AUDIO_STT_OPENAI_API_KEY` | str | `OPENAI_API_KEY` | **Yes** |
| `AUDIO_STT_SUPPORTED_CONTENT_TYPES` | str | `""` | **Yes** | Comma-separated content types |
| `AUDIO_STT_AZURE_API_KEY` | str | `""` | **Yes** |
| `AUDIO_STT_AZURE_REGION` | str | `""` | **Yes** |
| `AUDIO_STT_AZURE_LOCALES` | str | `""` | **Yes** |
| `AUDIO_STT_AZURE_BASE_URL` | str | `""` | **Yes** |
| `AUDIO_STT_AZURE_MAX_SPEAKERS` | str | `""` | **Yes** |
| `AUDIO_STT_MISTRAL_API_KEY` | str | `""` | **Yes** |
| `AUDIO_STT_MISTRAL_API_BASE_URL` | str | `https://api.mistral.ai/v1` | **Yes** |
| `AUDIO_STT_MISTRAL_USE_CHAT_COMPLETIONS` | bool | `False` | **Yes** |

### Text-to-Speech (TTS)

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `AUDIO_TTS_ENGINE` | str | `""` | **Yes** | TTS engine: `openai`, `elevenlabs`, `azure` |
| `AUDIO_TTS_MODEL` | str | `tts-1` | **Yes** | TTS model |
| `AUDIO_TTS_VOICE` | str | `alloy` | **Yes** | TTS voice |
| `AUDIO_TTS_SPLIT_ON` | str | `punctuation` | **Yes** | Split on: `punctuation`, `words`, `chars` |
| `AUDIO_TTS_API_KEY` | str | `""` | **Yes** | TTS API key |
| `ELEVENLABS_API_BASE_URL` | str | `https://api.elevenlabs.io` | No | ElevenLabs API URL |
| `AUDIO_TTS_OPENAI_API_BASE_URL` | str | `OPENAI_API_BASE_URL` | **Yes** |
| `AUDIO_TTS_OPENAI_API_KEY` | str | `OPENAI_API_KEY` | **Yes** |
| `AUDIO_TTS_OPENAI_PARAMS` | str | `{}` | **Yes** |
| `AUDIO_TTS_AZURE_SPEECH_REGION` | str | `""` | **Yes** |
| `AUDIO_TTS_AZURE_SPEECH_BASE_URL` | str | `""` | **Yes** |
| `AUDIO_TTS_AZURE_SPEECH_OUTPUT_FORMAT` | str | `audio-24khz-160kbitrate-mono-mp3` | **Yes** |

---

## User Permissions

### Workspace Permissions

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `USER_PERMISSIONS_WORKSPACE_MODELS_ACCESS` | bool | `False` | Access to workspace models |
| `USER_PERMISSIONS_WORKSPACE_KNOWLEDGE_ACCESS` | bool | `False` | Access to knowledge base |
| `USER_PERMISSIONS_WORKSPACE_PROMPTS_ACCESS` | bool | `False` | Access to prompts |
| `USER_PERMISSIONS_WORKSPACE_TOOLS_ACCESS` | bool | `False` | Access to tools |
| `USER_PERMISSIONS_WORKSPACE_SKILLS_ACCESS` | bool | `False` | Access to skills |
| `USER_PERMISSIONS_WORKSPACE_MODELS_IMPORT` | bool | `False` | Import models |
| `USER_PERMISSIONS_WORKSPACE_MODELS_EXPORT` | bool | `False` | Export models |
| `USER_PERMISSIONS_WORKSPACE_PROMPTS_IMPORT` | bool | `False` | Import prompts |
| `USER_PERMISSIONS_WORKSPACE_PROMPTS_EXPORT` | bool | `False` | Export prompts |
| `USER_PERMISSIONS_WORKSPACE_TOOLS_IMPORT` | bool | `False` | Import tools |
| `USER_PERMISSIONS_WORKSPACE_TOOLS_EXPORT` | bool | `False` | Export tools |

### Sharing Permissions

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `USER_PERMISSIONS_WORKSPACE_MODELS_ALLOW_SHARING` | bool | `False` | Share models |
| `USER_PERMISSIONS_WORKSPACE_MODELS_ALLOW_PUBLIC_SHARING` | bool | `False` | Publicly share models |
| `USER_PERMISSIONS_WORKSPACE_KNOWLEDGE_ALLOW_SHARING` | bool | `False` | Share knowledge |
| `USER_PERMISSIONS_WORKSPACE_KNOWLEDGE_ALLOW_PUBLIC_SHARING` | bool | `False` | Publicly share knowledge |
| `USER_PERMISSIONS_WORKSPACE_PROMPTS_ALLOW_SHARING` | bool | `False` | Share prompts |
| `USER_PERMISSIONS_WORKSPACE_PROMPTS_ALLOW_PUBLIC_SHARING` | bool | `False` | Publicly share prompts |
| `USER_PERMISSIONS_WORKSPACE_TOOLS_ALLOW_SHARING` | bool | `False` | Share tools |
| `USER_PERMISSIONS_WORKSPACE_TOOLS_ALLOW_PUBLIC_SHARING` | bool | `False` | Publicly share tools |
| `USER_PERMISSIONS_WORKSPACE_SKILLS_ALLOW_SHARING` | bool | `False` | Share skills |
| `USER_PERMISSIONS_WORKSPACE_SKILLS_ALLOW_PUBLIC_SHARING` | bool | `False` | Publicly share skills |
| `USER_PERMISSIONS_NOTES_ALLOW_SHARING` | bool | `False` | Share notes |
| `USER_PERMISSIONS_NOTES_ALLOW_PUBLIC_SHARING` | bool | `False` | Publicly share notes |

### Chat Permissions

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `USER_PERMISSIONS_CHAT_CONTROLS` | bool | `True` | Show chat controls |
| `USER_PERMISSIONS_CHAT_VALVES` | bool | `True` | Allow valves adjustment |
| `USER_PERMISSIONS_CHAT_SYSTEM_PROMPT` | bool | `True` | Allow system prompt editing |
| `USER_PERMISSIONS_CHAT_PARAMS` | bool | `True` | Allow parameter adjustment |
| `USER_PERMISSIONS_CHAT_FILE_UPLOAD` | bool | `True` | Allow file upload |
| `USER_PERMISSIONS_CHAT_DELETE` | bool | `True` | Allow chat deletion |
| `USER_PERMISSIONS_CHAT_DELETE_MESSAGE` | bool | `True` | Allow message deletion |
| `USER_PERMISSIONS_CHAT_CONTINUE_RESPONSE` | bool | `True` | Allow continue response |
| `USER_PERMISSIONS_CHAT_REGENERATE_RESPONSE` | bool | `True` | Allow regenerate |
| `USER_PERMISSIONS_CHAT_RATE_RESPONSE` | bool | `True` | Allow rating responses |
| `USER_PERMISSIONS_CHAT_EDIT` | bool | `True` | Allow message editing |
| `USER_PERMISSIONS_CHAT_SHARE` | bool | `True` | Allow chat sharing |
| `USER_PERMISSIONS_CHAT_EXPORT` | bool | `True` | Allow chat export |
| `USER_PERMISSIONS_CHAT_STT` | bool | `True` | Allow speech-to-text |
| `USER_PERMISSIONS_CHAT_TTS` | bool | `True` | Allow text-to-speech |
| `USER_PERMISSIONS_CHAT_CALL` | bool | `True` | Allow voice calls |
| `USER_PERMISSIONS_CHAT_MULTIPLE_MODELS` | bool | `True` | Allow multiple models |
| `USER_PERMISSIONS_CHAT_TEMPORARY` | bool | `True` | Allow temporary chats |
| `USER_PERMISSIONS_CHAT_TEMPORARY_ENFORCED` | bool | `False` | Force temporary chats |

### Feature Permissions

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `USER_PERMISSIONS_FEATURES_DIRECT_TOOL_SERVERS` | bool | `False` | Allow direct tool servers |
| `USER_PERMISSIONS_FEATURES_WEB_SEARCH` | bool | `True` | Allow web search |
| `USER_PERMISSIONS_FEATURES_IMAGE_GENERATION` | bool | `True` | Allow image generation |
| `USER_PERMISSIONS_FEATURES_CODE_INTERPRETER` | bool | `True` | Allow code interpreter |
| `USER_PERMISSIONS_FEATURES_FOLDERS` | bool | `True` | Allow folders |
| `USER_PERMISSIONS_FEATURES_NOTES` | bool | `True` | Allow notes |
| `USER_PERMISSIONS_FEATURES_CHANNELS` | bool | `True` | Allow channels |
| `USER_PERMISSIONS_FEATURES_API_KEYS` | bool | `False` | Allow API keys |
| `USER_PERMISSIONS_FEATURES_MEMORIES` | bool | `True` | Allow memories |
| `USER_PERMISSIONS_SETTINGS_INTERFACE` | bool | `True` | Allow interface settings |

---

## Task & Prompt Templates

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `TASK_MODEL` | str | `""` | **Yes** | Model for internal tasks (Ollama) |
| `TASK_MODEL_EXTERNAL` | str | `""` | **Yes** | Model for internal tasks (external APIs) |
| `ENABLE_TITLE_GENERATION` | bool | `True` | **Yes** | Auto-generate chat titles |
| `TITLE_GENERATION_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Custom title generation prompt |
| `ENABLE_TAGS_GENERATION` | bool | `True` | **Yes** | Auto-generate tags |
| `TAGS_GENERATION_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Custom tags prompt |
| `ENABLE_FOLLOW_UP_GENERATION` | bool | `True` | **Yes** | Generate follow-up questions |
| `FOLLOW_UP_GENERATION_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Custom follow-up prompt |
| `ENABLE_SEARCH_QUERY_GENERATION` | bool | `True` | **Yes** | Generate search queries |
| `ENABLE_RETRIEVAL_QUERY_GENERATION` | bool | `True` | **Yes** | Generate retrieval queries |
| `QUERY_GENERATION_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Custom query generation prompt |
| `ENABLE_AUTOCOMPLETE_GENERATION` | bool | `False` | **Yes** | Enable autocompletion |
| `AUTOCOMPLETE_GENERATION_INPUT_MAX_LENGTH` | int | `-1` | **Yes** | Max input length for autocomplete |
| `AUTOCOMPLETE_GENERATION_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Custom autocomplete prompt |
| `IMAGE_PROMPT_GENERATION_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Custom image prompt generation |
| `VOICE_MODE_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Voice mode system prompt |
| `TOOLS_FUNCTION_CALLING_PROMPT_TEMPLATE` | str | *(default)* | **Yes** | Tool calling prompt |

---

## Code Execution

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `ENABLE_CODE_EXECUTION` | bool | `True` | **Yes** | Enable code execution |
| `CODE_EXECUTION_ENGINE` | str | `pyodide` | **Yes** | Engine: `pyodide`, `jupyter` |
| `CODE_EXECUTION_JUPYTER_URL` | str | `""` | **Yes** | Jupyter server URL |
| `CODE_EXECUTION_JUPYTER_AUTH` | str | `""` | **Yes** | Auth method: `token`, `password` |
| `CODE_EXECUTION_JUPYTER_AUTH_TOKEN` | str | `""` | **Yes** | Jupyter auth token |
| `CODE_EXECUTION_JUPYTER_AUTH_PASSWORD` | str | `""` | **Yes** | Jupyter auth password |
| `CODE_EXECUTION_JUPYTER_TIMEOUT` | int | `60` | **Yes** | Execution timeout (seconds) |
| `ENABLE_CODE_INTERPRETER` | bool | `True` | **Yes** | Enable code interpreter |
| `CODE_INTERPRETER_ENGINE` | str | `pyodide` | **Yes** | Interpreter engine |
| `CODE_INTERPRETER_PROMPT_TEMPLATE` | str | `""` | **Yes** | Custom interpreter prompt |
| `CODE_INTERPRETER_JUPYTER_URL` | str | `CODE_EXECUTION_JUPYTER_URL` | **Yes** | Jupyter URL |
| `CODE_INTERPRETER_JUPYTER_AUTH` | str | `CODE_EXECUTION_JUPYTER_AUTH` | **Yes** | Auth method |
| `CODE_INTERPRETER_JUPYTER_AUTH_TOKEN` | str | `CODE_EXECUTION_JUPYTER_AUTH_TOKEN` | **Yes** | Auth token |
| `CODE_INTERPRETER_JUPYTER_AUTH_PASSWORD` | str | `CODE_EXECUTION_JUPYTER_AUTH_PASSWORD` | **Yes** | Auth password |
| `CODE_INTERPRETER_JUPYTER_TIMEOUT` | int | `60` | **Yes** | Timeout |
| `CODE_INTERPRETER_BLOCKED_MODULES` | str | `""` | No | Comma-separated blocked Python modules |

---

## Storage Providers

### Local

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `DATA_DIR` | str | `./data` | Base data directory |
| `UPLOAD_DIR` | str | `{DATA_DIR}/uploads` | Upload directory |
| `CACHE_DIR` | str | `{DATA_DIR}/cache` | Cache directory |
| `STATIC_DIR` | str | `./static` | Static files directory |
| `FRONTEND_BUILD_DIR` | str | `../build` | Frontend build directory |

### S3 Storage

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `STORAGE_PROVIDER` | str | `local` | Provider: `local`, `s3`, `gcs`, `azure` |
| `S3_ACCESS_KEY_ID` | str | `None` | AWS access key |
| `S3_SECRET_ACCESS_KEY` | str | `None` | AWS secret key |
| `S3_REGION_NAME` | str | `None` | AWS region |
| `S3_BUCKET_NAME` | str | `None` | S3 bucket name |
| `S3_KEY_PREFIX` | str | `None` | Key prefix for objects |
| `S3_ENDPOINT_URL` | str | `None` | Custom endpoint (MinIO, etc.) |
| `S3_USE_ACCELERATE_ENDPOINT` | bool | `False` | Use S3 Transfer Acceleration |
| `S3_ADDRESSING_STYLE` | str | `None` | Addressing style: `path`, `virtual` |
| `S3_ENABLE_TAGGING` | bool | `False` | Enable object tagging |

### Google Cloud Storage

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `GCS_BUCKET_NAME` | str | `None` | GCS bucket name |
| `GOOGLE_APPLICATION_CREDENTIALS_JSON` | str | `None` | Service account JSON |

### Azure Storage

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `AZURE_STORAGE_ENDPOINT` | str | `None` | Azure Storage endpoint |
| `AZURE_STORAGE_CONTAINER_NAME` | str | `None` | Container name |
| `AZURE_STORAGE_KEY` | str | `None` | Storage account key |

---

## Logging & Monitoring

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `GLOBAL_LOG_LEVEL` | str | `INFO` | Log level: `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL` |
| `ENABLE_AUDIT_STDOUT` | bool | `False` | Output audit logs to stdout |
| `ENABLE_AUDIT_LOGS_FILE` | bool | `True` | Write audit logs to file |
| `AUDIT_LOGS_FILE_PATH` | str | `{DATA_DIR}/audit.log` | Audit log file path |
| `AUDIT_LOG_FILE_ROTATION_SIZE` | str | `10MB` | Rotation size (e.g., `10MB`, `1GB`) |
| `AUDIT_UVICORN_LOGGER_NAMES` | str | `uvicorn.access` | Loggers to capture |
| `AUDIT_LOG_LEVEL` | str | `NONE` | Audit detail: `NONE`, `METADATA`, `REQUEST`, `REQUEST_RESPONSE` |
| `MAX_BODY_LOG_SIZE` | int | `2048` | Max body size to log (bytes) |
| `AUDIT_EXCLUDED_PATHS` | str | `/chats,/chat,/folders` | Paths to exclude from audit |

---

## Performance Tuning

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `THREAD_POOL_SIZE` | int | `0` (40) | FastAPI thread pool. **Increase for large instances** |
| `UVICORN_WORKERS` | int | `1` | Number of Uvicorn workers |
| `ENABLE_ASYNC_EMBEDDING` | bool | `True` | Async embedding generation |
| `RAG_EMBEDDING_BATCH_SIZE` | int | `1` | Embedding batch size |
| `CHAT_RESPONSE_STREAM_DELTA_CHUNK_SIZE` | int | `1` | Min tokens to batch before sending |
| `CHAT_STREAM_RESPONSE_CHUNK_MAX_BUFFER_SIZE` | int | `None` | Max chunk buffer size (bytes) |
| `AIOHTTP_CLIENT_TIMEOUT` | int | `None` (300) | HTTP client timeout |
| `AIOHTTP_CLIENT_TIMEOUT_MODEL_LIST` | int | `10` | Timeout for model list fetching |
| `AIOHTTP_CLIENT_TIMEOUT_TOOL_SERVER_DATA` | int | `10` | Timeout for tool server |
| `DATABASE_POOL_SIZE` | int | `None` | Database connection pool size |
| `DATABASE_POOL_MAX_OVERFLOW` | int | `0` | Max overflow connections |
| `DATABASE_POOL_TIMEOUT` | int | `30` | Pool timeout (seconds) |
| `DATABASE_POOL_RECYCLE` | int | `3600` | Connection recycle (seconds) |
| `ENABLE_COMPRESSION_MIDDLEWARE` | bool | `True` | Gzip compression |
| `MODELS_CACHE_TTL` | int | `1` | Model list cache TTL (seconds) |

---

## WebSocket Configuration

| Variable | Type | Default | Description |
|----------|------|---------|-------------|
| `ENABLE_WEBSOCKET_SUPPORT` | bool | `True` | Enable WebSocket support |
| `WEBSOCKET_MANAGER` | str | `""` | WebSocket manager implementation |
| `WEBSOCKET_REDIS_URL` | str | `REDIS_URL` | Redis URL for WebSocket pub/sub |
| `WEBSOCKET_REDIS_CLUSTER` | bool | `REDIS_CLUSTER` | Use Redis Cluster |
| `WEBSOCKET_REDIS_OPTIONS` | str | `""` | JSON options for Redis connection |
| `WEBSOCKET_SENTINEL_HOSTS` | str | `""` | Sentinel hosts |
| `WEBSOCKET_SENTINEL_PORT` | str | `26379` | Sentinel port |
| `WEBSOCKET_SERVER_LOGGING` | bool | `False` | Enable WebSocket server logging |
| `WEBSOCKET_SERVER_ENGINEIO_LOGGING` | bool | `False` | Enable Engine.IO logging |
| `WEBSOCKET_SERVER_PING_TIMEOUT` | int | `20` | Ping timeout (seconds) |
| `WEBSOCKET_SERVER_PING_INTERVAL` | int | `25` | Ping interval (seconds) |
| `WEBSOCKET_REDIS_LOCK_TIMEOUT` | int | `60` | Lock timeout (seconds) |

---

## Tool Servers

| Variable | Type | Default | Persistent | Description |
|----------|------|---------|------------|-------------|
| `TOOL_SERVER_CONNECTIONS` | str | `[]` | **Yes** | JSON array of tool server connections |

---

## Deprecated Variables

| Deprecated | Replacement | Notes |
|------------|-------------|-------|
| `OLLAMA_API_BASE_URL` | `OLLAMA_BASE_URL` | Use OLLAMA_BASE_URL instead |
| `ENABLE_API_KEY` | `ENABLE_API_KEYS` | Plural form |
| `ENABLE_API_KEY_ENDPOINT_RESTRICTIONS` | `ENABLE_API_KEYS_ENDPOINT_RESTRICTIONS` | Plural form |
| `API_KEY_ALLOWED_ENDPOINTS` | `API_KEYS_ALLOWED_ENDPOINTS` | Plural form |
| `WEBUI_JWT_SECRET_KEY` | `WEBUI_SECRET_KEY` | Consolidated naming |
| `ENABLE_ADMIN_WORKSPACE_CONTENT_ACCESS` | `BYPASS_ADMIN_ACCESS_CONTROL` | More descriptive name |
| `OAUTH_GROUP_CLAIM` | `OAUTH_GROUPS_CLAIM` | Plural form |

---

## Security Recommendations

### Production Checklist

1. **Change Default Secrets:**
   ```bash
   WEBUI_SECRET_KEY=<strong-random-key>
   ```

2. **Enable Authentication:**
   ```bash
   WEBUI_AUTH=True
   ENABLE_PASSWORD_AUTH=True  # Keep until OAuth is verified
   ```

3. **Set Correct URL:**
   ```bash
   WEBUI_URL=https://your-domain.com  # Must be set BEFORE OAuth
   ```

4. **CORS Configuration:**
   ```bash
   CORS_ALLOW_ORIGIN=https://your-domain.com  # NOT *
   ```

5. **Disable Dangerous Features:**
   ```bash
   JWT_EXPIRES_IN=4w  # Never use -1 in production
   ENABLE_RAG_LOCAL_WEB_FETCH=False
   ```

6. **Database Security:**
   ```bash
   DATABASE_URL=postgresql://user:pass@host/db  # Use Postgres for production
   ```

7. **OAuth Security:**
   ```bash
   ENABLE_PASSWORD_AUTH=False  # Only after OAuth is working
   OAUTH_ALLOWED_DOMAINS=yourcompany.com
   ```

---

## Multi-Instance Deployment

For scaling across multiple instances, these **must be identical** on all nodes:

| Variable | Why |
|----------|-----|
| `WEBUI_SECRET_KEY` | JWT tokens signed by one node must be valid on others |
| `REDIS_URL` | Shared state for WebSockets and config |
| `DATABASE_URL` | Shared database |
| `VECTOR_DB` + Vector DB config | Shared vector store |

### Required for Multi-Instance

```bash
# Redis is required for multi-instance
REDIS_URL=redis://redis:6379/0

# WebSocket pub/sub
WEBSOCKET_REDIS_URL=redis://redis:6379/0

# Database (PostgreSQL recommended)
DATABASE_URL=postgresql://user:pass@postgres:5432/openwebui

# Same secret on all nodes
WEBUI_SECRET_KEY=<same-random-key>
```

---

---

**OpenWebUI Version:** v0.8.3
**Last Updated:** 2026-02-18
