# Security Documentation

Last updated: 2026-02-18 (UTC)

## Security Policy

### Reporting Vulnerabilities

If you discover a security vulnerability in this deployment or its dependencies, please report it responsibly:

- **Email**: security@example.com
- **Response Time**: We aim to acknowledge reports within 48 hours
- **Disclosure Policy**: Please allow reasonable time for remediation before public disclosure

### What to Include

When reporting, please provide:
- Description of the vulnerability
- Steps to reproduce
- Potential impact assessment
- Suggested fix (if available)

## Supported Versions

| Component | Current Version | Support Status |
|-----------|-----------------|----------------|
| OpenWebUI | v0.8.3 | Active |
| CLIProxyAPI | v6.8.18 | Active |
| pgvector (PostgreSQL 16) | pg16 | Active |
| Jupyter | Latest | Active |

## Known Vulnerabilities

### Current CVE Status by Image

| Image | Version | CRITICAL | HIGH | MEDIUM | LOW | Status |
|-------|---------|----------|------|--------|-----|--------|
| OpenWebUI | v0.8.3 | 5 | 178 | - | - | Under Review |
| CLIProxyAPI | v6.8.18 | 0 | - | - | - | Updated (was v6.2.38 with 2 CRITICAL) |
| pgvector | pg16 | 2 | 6 | - | - | Under Review |
| Jupyter | Latest | - | - | - | - | Updated |

### Detailed Findings

#### OpenWebUI v0.8.3
- **5 CRITICAL** CVEs pending remediation
- **178 HIGH** CVEs pending triage
- Action Required: Evaluate upstream patches or image rebuild

#### CLIProxyAPI v6.8.18
- **Status**: Resolved
- Previous version (v6.2.38) had 2 CRITICAL CVEs
- Updated to v6.8.18 with vulnerabilities addressed

#### pgvector:pg16
- **2 CRITICAL** CVEs under investigation
- **6 HIGH** CVEs under investigation
- Action Required: Monitor PostgreSQL security advisories

#### Jupyter
- **Status**: Updated to recent version
- No pending critical vulnerabilities

## Acceptance Criteria

### Response Timeframes

| Severity | Maximum Response Time | Action Required |
|----------|----------------------|-----------------|
| CRITICAL | 7 days | Immediate patching or mitigation required |
| HIGH | 30 days | Patching schedule required |
| MEDIUM | 90 days | Include in regular update cycle |
| LOW | Next maintenance window | Address as convenient |

### Escalation Path

1. **CRITICAL**: Immediate notification to infrastructure team
2. **HIGH**: Security review within 7 days
3. **MEDIUM/LOW**: Document and schedule

## Update History

| Date | Component | Action | Result |
|------|-----------|--------|--------|
| 2026-02-18 | CLIProxyAPI | v6.2.38 -> v6.8.18 | Resolved 2 CRITICAL CVEs |
| 2026-02-18 | Jupyter | Updated to latest | Security hardening |
| 2026-02-18 | SECURITY.md | Created | Baseline documentation |

---

## Vulnerability Scanning

To scan current images for vulnerabilities:

```bash
# Using Trivy (recommended)
trivy image ghcr.io/open-webui/open-webui:0.8.3

# Using Docker Scout
docker scout cves ghcr.io/open-webui/open-webui:0.8.3
```

## References

- [OpenWebUI Security Advisories](https://github.com/open-webui/open-webui/security)
- [PostgreSQL Security](https://www.postgresql.org/support/security/)
- [NIST CVE Database](https://nvd.nist.gov/)
