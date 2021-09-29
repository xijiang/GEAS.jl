## Notes to the codes in `tst/`
These codes will be deleted in the final release, except:

1. Example codes for parameter setup.
2. Example codes for `slurm` tasks.


## Log
### 2021-09-23
Developing base moved from a3970x to 3900x.
Saga is much more powerful with so many cores and memory.

- Problems
  - in function `create_storage`, `:tbv`, and `p7e` were not calculated in generation 0.
