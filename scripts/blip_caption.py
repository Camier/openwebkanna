#!/usr/bin/env python3
"""Generate an image caption using Salesforce BLIP‑Base.

The script loads the model "salesforce/blip-image-captioning-base" via the
transformers library, runs inference on a supplied image file, and prints the
caption to STDOUT. It is deliberately minimal – the surrounding Bash wrapper
handles iteration, error handling and persistence.

Dependencies (install in the container where you run the script):
    pip install torch torchvision transformers pillow
"""

import sys
from pathlib import Path

import torch
from PIL import Image
from transformers import BlipProcessor, BlipForConditionalGeneration


def load_model():
    """Load the BLIP‑Base processor and model.

    Returns a tuple (processor, model) ready for inference.
    """
    processor = BlipProcessor.from_pretrained("salesforce/blip-image-captioning-base")
    model = BlipForConditionalGeneration.from_pretrained(
        "salesforce/blip-image-captioning-base"
    )
    # Use mixed‑precision if a GPU is available
    if torch.cuda.is_available():
        model.to("cuda")
        torch.set_float32_matmul_precision("high")
    return processor, model


def caption_image(image_path: Path, processor, model) -> str:
    """Generate a caption for *image_path*.

    The function returns the generated text string. If inference fails it
    returns an empty string.
    """
    try:
        raw_image = Image.open(image_path).convert("RGB")
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open image {image_path}: {e}\n")
        return ""

    inputs = processor(images=raw_image, return_tensors="pt")
    if torch.cuda.is_available():
        inputs = {k: v.to("cuda") for k, v in inputs.items()}
    try:
        out = model.generate(**inputs, max_length=64, num_beams=4)
        caption = processor.decode(out[0], skip_special_tokens=True)
        return caption
    except Exception as e:
        sys.stderr.write(f"[ERROR] Inference failed for {image_path}: {e}\n")
        return ""


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: blip_caption.py <image_path>\n")
        sys.exit(1)
    image_path = Path(sys.argv[1])
    if not image_path.is_file():
        sys.stderr.write(f"[ERROR] File not found: {image_path}\n")
        sys.exit(1)

    processor, model = load_model()
    caption = caption_image(image_path, processor, model)
    if caption:
        print(caption)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
