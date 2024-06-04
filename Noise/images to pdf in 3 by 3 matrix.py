from PIL import Image
import os

def resize_image(image, target_width, target_height):
    """
    Resize an image to the target width and height while maintaining aspect ratio.
    """
    return image.resize((target_width, target_height), Image.ANTIALIAS)

def merge_images_into_pdf(folder_path, output_pdf_path):
    image_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.jpg') or f.endswith('.png')])

    images = {'normalized_psd_plot_': [], 'psd_plot_slope_it_': [], 'signal_after_decimation_': []}
    for file in image_files:
        category = None
        for key in images.keys():
            if file.startswith(key):
                category = key
                break

        if category is not None:
            image_path = os.path.join(folder_path, file)
            images[category].append(resize_image(Image.open(image_path), int(4.3 * 100), int(2.5 * 100)))  # 100 DPI
        else:
            print(f"Ignoring file {file} because it doesn't match any category.")

    if not all(images.values()):
        print("Some categories do not have images.")
        return

    pdf_path = output_pdf_path
    num_images_per_page = 3

    # Create a list to store all pages
    pages = []

    # Iterate over the images and create pages
    for i in range(0, len(images['normalized_psd_plot_']), num_images_per_page):
        # Create a new page
        page = Image.new('RGB', (3 * int(4.3 * 100), 3 * int(2.5 * 100)))  # 100 DPI

        # Paste three images from each category onto the page
        for j in range(num_images_per_page):
            for k, category in enumerate(images.keys()):
                img_index = i + j
                img = images[category][img_index] if img_index < len(images[category]) else None
                if img:
                    col = k
                    row = j
                    x_offset = col * img.width
                    y_offset = row * img.height
                    page.paste(img, (x_offset, y_offset))

        pages.append(page)

    # Save all pages into a single PDF
    pages[0].save(pdf_path, "PDF", resolution=100.0, save_all=True, append_images=pages[1:], quality=95)

    print("PDF created successfully.")

if __name__ == "__main__":
    folder_path = r"C:\Users\sndkp\OneDrive\Desktop\plots new 0.0001 to 0.1 - Copy"
    output_pdf_path = "output1.pdf"
    merge_images_into_pdf(folder_path, output_pdf_path)
