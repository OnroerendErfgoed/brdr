import os

from examples import save_animated_gif

output_filename = "test.gif"
frame_duration_ms = 500

image_files = sorted(
    [f for f in os.listdir(".") if f.endswith(".png") or f.endswith(".jpg")]
)
save_animated_gif(image_files, output_filename, frame_duration_ms)
