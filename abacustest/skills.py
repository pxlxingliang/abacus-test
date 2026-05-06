import os


def show_skills(param=None):
    main_dir = os.path.dirname(os.path.abspath(__file__))
    skills_path = os.path.join(main_dir, "skills")

    if not os.path.isdir(skills_path):
        print(f"ERROR: skills directory not found at {skills_path}")
        return

    print(f"abacustest skills path: {skills_path}")
    print("\nThis is the abacustest skills path. You can add it to the skills directory")
    print("of openclaw/opencode/clawcode and other AI agents. After that, your AI agent")
    print("will be able to correctly handle ABACUS tasks and use abacustest.")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="abacustest skills")
    param = parser.parse_args()
    show_skills(param)


if __name__ == "__main__":
    main()