ffmpeg -framerate 60 -i frame%04d.png output.mp4 -c:a copy -c:v libx264 -crf 18