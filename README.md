# Morpheus

Probing whether single-cell foundation models (Geneformer) implicitly encode known regulatory hierarchies.

See [`CLAUDE.md`](./CLAUDE.md) for the project overview, [`shared/contracts.md`](./shared/contracts.md) for cross-component data shapes, and [`shared/tickets-status.md`](./shared/tickets-status.md) for the backlog and progress log.

## Quick start

```bash
uv sync --extra dev
uv run pytest
```

Then follow the staged commands in `CLAUDE.md` (data → embed → probe → analysis).
