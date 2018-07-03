To run the clustering analysis on O2, be sure to have X11 forwarding working if you want to visualize any of the images. To do this, you may need to have XQuartz running on your local machine and log onto O2 with the terminal:

```bash
ssh -XY username@o2.hms.harvard.edu
```

Then start an interactive session with extra memory and x11:

```bash
srun --pty -p interactive -t 0-12:00 --x11 --mem 64G /bin/bash
```

After starting the interactive session, load the necessary R modules:

```bash
module load gcc/6.2.0 R/3.4.1
```

