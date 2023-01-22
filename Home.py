import streamlit as st
from streamlit_option_menu import option_menu

import base64


def body_bg(side_bg):
    '''
    A function to unpack an image from url and set as bg.
    Returns
    -------
    The background.
    '''
    side_bg_ext = 'png'

    st.markdown(
        f"""
         <style>
         .stApp {{
            background: url(data:image/{side_bg_ext};base64,{base64.b64encode(open(side_bg, "rb").read()).decode()});
            background-size: cover
         }}
         </style>
         """,
        unsafe_allow_html=True
    )


def sidebar_bg(side_bg):

   side_bg_ext = 'png'

   st.markdown(
      f"""
      <style>
      [data-testid="stBody"] > div:first-child {{
          background: url(data:image/{side_bg_ext};base64,{base64.b64encode(open(side_bg, "rb").read()).decode()});
      }}
      </style>
      """,
      unsafe_allow_html=True,
      )

body_bg('./images/background2.jpeg')
sidebar_bg('./images/background2.jpeg')


selected = option_menu(
    menu_title=None,
    options = ["Home", "Technology", "Team"],
    # icons = []
    default_index = 0,
    orientation = "horizontal",
)

if selected == "Home":
    st.title(f"You have selected {selected}")
if selected == "Technology":
    st.title(f"You have selected {selected}")
if selected == "Team":
    st.title(f"You have selected {selected}")

if selected == "Demo":
    st.title(f"You have selected {selected}")