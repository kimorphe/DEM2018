ffmpeg \
-i movie.mp4 \
-i Rstat1/movie.mp4 \
-i Rstat2/movie.mp4 \
-i Rstat3/movie.mp4 \
-i Rstat4/movie.mp4 \
-i Rstat5/movie.mp4 \
-i Rstat6/movie.mp4 \
-i Rstat7/movie.mp4 \
-i Rstat8/movie.mp4 \
 -strict -2 -filter_complex 'concat=n=9:v=1' mvs.mp4
