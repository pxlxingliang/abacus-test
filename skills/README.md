# ABACUS Test Skills

This repository contains a collection of **AgentSkills** designed for ABACUS (Atomic-orbital Based Ab-initio Computation at University of Science and technology of China) workflows. These skills enable AI agents (OpenClaw, ClawCode, OpenCode, etc.) to interact with ABACUS calculations programmatically.

## 📦 What Are Skills?

**AgentSkills** are modular, reusable instruction sets that extend AI agent capabilities. Each skill:

- Defines **when** to use it (trigger conditions)
- Specifies **how** to execute tasks (commands, parameters, workflows)
- Provides **examples** and **limitations**
- Can be loaded by compatible AI agent frameworks

Think of skills as "plugins" that give your AI agent domain-specific expertise in ABACUS computational materials science workflows.


## 🔌 Integration Guide

### For OpenClaw

1. **Copy skills to your workspace:**
   ```bash
   cp -r skills/* /root/.openclaw/workspace/skills/
   ```

2. **Register in OpenClaw config** (if required by your version):
   ```json
   {
     "skills": {
       "paths": ["/root/.openclaw/workspace/skills"]
     }
   }
   ```

3. **Restart OpenClaw gateway:**
   ```bash
   openclaw gateway restart
   ```

4. **Verify skills are loaded:**
   - Ask your agent: "What skills do you have?"
   - Or check: `openclaw skills list`

### For ClawCode / OpenCode

1. **Add skills path to agent configuration:**
   ```json
   {
     "agent": {
       "skills": {
         "enabled": true,
         "paths": ["PATH-TO-SKILLS/skills"]
       }
     }
   }
   ```

2. **Or symlink skills:**
   ```bash
   ln -s /ABSOLUTE-PATH-TO-SKILLS/skills/* ~/.clawcode/skills/
   ```

3. **Reload agent context** (method varies by platform)

### General Integration (Any Agent Framework)

Skills follow the **AgentSkills specification**:

- Each skill is a directory containing `SKILL.md`
- `SKILL.md` has YAML frontmatter + markdown content
- Frontmatter includes: `name`, `description`, `metadata`

**Minimum requirements:**
```yaml
---
name: skill-name
description: "What this skill does"
metadata: { "openclaw": { "emoji": "🔧" } }
---
```

## 🚀 Quick Start

Once skills are loaded, you can use them naturally:

**Examples:**

> "Convert this POSCAR to ABACUS STRU format"
> → Triggers `abacustest-prepare`

> "Extract the total energy from this calculation"
> → Triggers `abacustest-extract-dft-results`

> "Convert my ABACUS inputs to VASP format"
> → Triggers `abacustest-abacus2VaspQeCp2k`

The agent will automatically select the appropriate skill based on your request.

## 📋 Requirements

- **Python package**: `abacustest` (install via `pip install abacustest`)
- **Compatible agent framework**: OpenClaw, ClawCode, OpenCode, or any AgentSkills-compatible system

