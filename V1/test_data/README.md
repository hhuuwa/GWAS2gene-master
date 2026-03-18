# V1 Test Data

This folder provides lightweight, shareable demo data and expected outputs.

## Included files

- `RAP_HZ_Awn_length_demo.FarmCPUpeak_info`
  - Small subset extracted from the full HZ demo intermediate table.
- `RAP_Step2_HZ_Awn_length_demo.FarmCPUpeak_info`
  - Step2 result generated from the subset above.

## Why this is useful

- Readers can quickly understand expected file format.
- You can run a fast sanity check without full data volume.

## Fast check command

```powershell
<python> V1\code\RAP_Stpe2_select_FarmCPUpeak_info.py HZ_Awn_length_demo Awn_length
```
