# WRF-GRK Agent Notes

Before doing WRF-GRK work, read these files in order:

1. `docs/agent_memory.md` - current operational memory for agent sessions.
2. `docs/simulation.md` - authoritative scientific/run history and bug log.
3. `/home/igrk/.agents/skills/wrf-ghg/SKILL.md` - workflow rules and helper scripts.
4. `/home/igrk/.agents/skills/wrf-ghg/references/segment_workflow.md` - segment launch checklist.

Important defaults:

- Work from `/home/igrk/WRF-GRK`.
- Use `mpirun -np 16` for WRF/real on this host.
- Do not launch a segment unless disk free space is safely above the segment budget.
- Always verify segment inputs with:

```bash
.venv/bin/python /home/igrk/.agents/skills/wrf-ghg/scripts/check_seg_inputs.py <YYYY-MM-DD start> <YYYY-MM-DD end>
```

- Use `Times` inside `wrfout` files as the time truth. Restart segment filenames begin at `_01:00:00`.
- Treat `docs/agent_memory.md` as the current live operational state. `docs/simulation.md` remains the science/bug-history reference, but some progress-log entries may describe older run attempts.
