# Stage I/O Contract

## xTB batch input

```csv
idx,name,xyz_path
```

## xTB batch output

```csv
idx,name,xyz_path,status,detail,total_energy_eh,homo_lumo_gap_ev,normal_termination
```

## TDDFT-xTB filter output

```csv
idx,name,xyz_path,status,detail,s1_ev,lambda_nm,target_nm,delta_nm,pass_blue_window
```

## Required logs

- `xtb_state.json`
- `xtb_progress.csv`
- `remote_xtb.log` (for remote runs)
