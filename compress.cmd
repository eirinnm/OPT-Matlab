FOR %%A IN (%*) DO (
"%~dp0\ffmpeg.exe" -i %%~A -acodec copy -vcodec ffv1 -level 3 -threads 4 -coder 1 -context 1 %%~dpnA_compressed.avi
)