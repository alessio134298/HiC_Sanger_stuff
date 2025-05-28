#!/usr/bin/env python3

import os
import sys
import json
import re  # Import regular expressions for parsing

from google import genai
from google.genai import types

import pymupdf4llm


FILE_DIR = 'test_files'
GLOBAL_DICT = {}

# Set up API key
print("Obtaining API key from environment variables.")
api_key = os.getenv("GOOGLE_API_KEY")
if api_key is None:
    print("❌ No API key found. Please set env variable 'GOOGLE_API_KEY'")
    sys.exit(1)
client = genai.Client(api_key=api_key)
print("✅ API Configured Successfully!")


# Set up model
MODEL_ID = 'gemini-2.5-flash-preview-04-17'
print(f"\nUsing model: '{MODEL_ID}'")

system_instructions = """You are an expert research assistant.
Your task is to classify each paper in the directory according to up to five of the
pre-specified categories. Please use the available tools to:
    1) find the list of PDFs in the directory,
    2) read the PDFs in the directories as markdown,
    3) classify each PDF to up to five categories.
    
Then, print the result."""


def list_files_in_dir() -> str:  # Return string for easier handling
    """Lists files in the file directory."""
    try:
        files = os.listdir(f'./{FILE_DIR}')
        if not files:
            return f"The '{FILE_DIR}' directory is empty."
        return "Files found: " + ", ".join(files)
    except FileNotFoundError:
        return f"Error: '{FILE_DIR}' directory not found."
    except Exception as e:
        return f"Error listing {FILE_DIR}: {e}"


def get_pdf_markdown(filename: str) -> str:
    """Reads a PDF file in the files directory and returns markdown.

    Args:
      filename: Filename of PDF without folder structure. E.g. paper1.pdf
    """
    pdf_path = os.path.join(f'{FILE_DIR}', filename)
    if not os.path.exists(pdf_path):
        return f"Error: PDF file not found at '{pdf_path}'"
    else:
        try:
            print(f"--- Processing PDF: {pdf_path} ---")
            md_text = pymupdf4llm.to_markdown(
                pdf_path, write_images=False, embed_images=False
            )
            print(f"--- Finished processing PDF: {pdf_path} ---")
            if not md_text or md_text.isspace():
                return f"Successfully processed PDF '{filename}', but it appears to contain no extractable text."
            return md_text
        except Exception as e:
            if "cannot open broken document" in str(e).lower():
                return f"Error processing PDF '{filename}': The file appears to be corrupted or is not a valid PDF."
            else:
                return f"An error occurred during PDF processing for '{filename}': {e}"


def store_paper_categories(categories: list[str], key: str) -> None:
    """
    Given a list of up to five categories (list[str]), stores these in a global
    dict under key 'key'. I would make the key the filename of the PDF.
    """
    GLOBAL_DICT[key] = categories
    return None


def get_all_categories() -> dict:
    return GLOBAL_DICT


# for gemini 2.5
# gen_config = types.GenerateContentConfig(
#     tools=[list_files_in_dir,
#            get_pdf_markdown,
#            store_paper_categories,
#            get_all_categories],
#     temperature=.5,
#     system_instruction=system_instructions,
#     thinking_config=types.ThinkingConfig(
#         includeThoughts=True, thinking_budget=0)
# )

gen_config = types.GenerateContentConfig(
    tools=[list_files_in_dir,
           get_pdf_markdown,
           store_paper_categories,
           get_all_categories],
    temperature=.5,
    system_instruction=system_instructions
)


chat = client.chats.create(
    model=MODEL_ID,
    config=gen_config
)

prompt = """
Please perform the following actions using the available tools:
1. Find all PDFs of papers in the directory.
2. Read each paper and get the markdown.
3. Classify the papers into pre-set categories, store the result in the dictionary. (coming up)
4. Print the result.

If there are many papers, please classify each paper immediately after reading.

Your categories are these:
[Bioinformatics, Genomics, Psychology, Immunology, Racing, Playmobile, Barbecue]
"""

print("\n--- Sending Prompt to LLM ---")
print(prompt)
print("-----------------------------")

response = chat.send_message(prompt)


print("\n--- Full Chat History ---")
for content in chat.get_history():
    print("### " + content.role.upper() + ":")  # Use upper for clarity
    for part in content.parts:
        if part.text:
            print(part.text)
        if fc := part.function_call:  # Use assignment expression for brevity
            print(f"Tool Call: {fc.name}({fc.args})")
        if fr := part.function_response:  # Use assignment expression
            # Shorten output for readability if response is very long (e.g. PDF markdown)
            response_content = str(fr.response)
            if len(response_content) > 500:
                response_content = response_content[:250] + \
                    "\n...\n" + response_content[-250:]
            print(f"Tool Response ({fr.name}): {response_content}")
            # print(f"Tool Response: {fr.name}(response={fr.response})") # Original more verbose format
    print("-" * 80)
