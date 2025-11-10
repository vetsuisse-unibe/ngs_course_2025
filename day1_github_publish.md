# 🚀 Publishing & Sharing Your Work

In this session, you’ll learn how to **prepare results**, **document your analysis**, **version your work with Git**, and **publish to GitHub Pages** so others can access your materials.

---

### 🧠 Learning Objectives

- Organize outputs and write clear documentation (README/Markdown).
- Use Git to version, tag, and sync your repository with GitHub.
- Publish a simple website with **GitHub Pages**.
- Share your work responsibly with **licenses**, **citations**, and **releases/DOIs**.

---

## 1️⃣ Repository Structure

- `hpc-exercises/`: Contains HPC cluster practice exercises
  - `scripts/`: Shell scripts for various HPC tasks
- [Additional directories will be added as we progress]

---

## 2️⃣ Exercise 1: GitHub Setup

### Create a New GitHub Repository
1. **Go to GitHub and Sign In:** Navigate to github.com and enter your username and password. If you are a first time user, sign up and login.
2. **Start a New Repository:**
  - Click the plus sign (+) icon located in the top-right corner of the page.
  - Select New repository from the dropdown menu.
3. **Choose a Good Name:**
  - In the "Repository name" field, enter a name that clearly describes the content of the project.
  - Good Examples: bioinformatics-exercises, Sequencing-course-exercises,etc
  - Avoid: Names like my-Project or test-repo
4. **Set Visibility:** Leave the repository set to Public.
5. **Skip the README:** Ensure the option "Add a README file" is unchecked. (You'll add files later.)
6. **Final Step:** Click the green Create repository button.
7. **Get the SSH Address:** On the next screen, find the area showing the repository address. Make sure the SSH tab is selected, and then copy the URL provided (it will start with git@github.com:...).

You just created your first scientific notebook on GitHub!. Feel free to store all your other code and projects here too - it's like having unlimited digital lab notebooks with automatic backup. **Generally, keep one GitHub repository per analysis project.**

### Configure Git on HPC Cluster (Bioinformatics server)
- Use the remote login extension on Visual studio code and login into the bioinformatics server (login8.hpc.binf.unibe.ch) with your chosen username and password.

#### Set Up Authentication (Important!)
To securely link your local git repository to GitHub, we use SSH authentication. Think of SSH as a digital keycard system:
- You create a unique key pair on your computer.

- You give GitHub the Public Key (like registering your keycard).

- Every time you try to push new scripts or large processed datasets to your repository, GitHub instantly checks your local computer's Private Key to verify your identity.

This method is more secure than using a password and makes updating GitHub much faster and easier for both scripts and high-volume data files.

```shell
# Generate SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"
```
-t ed25519: Specifies the type of key (Ed25519 is a modern, secure algorithm)

When you run this, it will:
- Ask where to save the key. Default location for eg. (/home/\<username\>/.ssh/id_ed25519)is good so just press _ENTER_
- Ask for a passphrase (optional but recommended for security)
A passphrase is like a password but typically longer and used to encrypt your SSH private key. If **you don't want to set it just type _ENTER_**. If you set please make sure you remember it or store it somewhere. **If you set a passphrase, you will be prompted to enter it every time you push to GitHub.**

```shell
# Display your public key (copy this output)
cat ~/.ssh/id_ed25519.pub
```

#### Linking Your Secure Key to GitHub
This process allows GitHub to recognize your computer instantly and securely when you try to upload files.

1. **Open GitHub Settings:**
  - Navigate to GitHub.com and make sure you are logged in.
  - Click your profile picture (top-right corner), and select Settings.

2. **Go to Key Management:**
  - In the left sidebar menu, find and click SSH and GPG keys.
  - Click the New SSH key button

3. **Upload Your Key:**
  - Give your key a descriptive Title (e.g., "My Laptop" or "HPC Access").
  - In the large "Key" field, paste the entire contents of your Public SSH Key. (cat ~/.ssh/id_ed25519.pub)
  - Click the green Add SSH key button to save it.
4. **Use the Right Address:**
  - When you connect your local project to GitHub, always use the SSH URL format (it starts with git@github.com:...) instead of the standard HTTPS URL.

---

## 3️⃣ Exercise 2: Connect Local git repository to Remote GitHub.

### Link Your HPC Repository
```shell
# Navigate to your existing repository
cd hpc-exercises

# Add remote repository (replace 'git@github.com:<username>/sequencencing-exercises.git' with your URL that copy it from your github repository page )
git remote add origin git@github.com:<username>/sequencencing-exercises.git

# Verify remote was added
git remote -v
```
The `git remote -v` command shows all the remote repositories connected to your local repository, along with their URLs. The -v stands for "verbose", showing both fetch and push URLs.

---

## 4️⃣ Exercise 3: Document Your Repository

### Create README.md
Create a new file called `README.md` in your repository root or project folder (hpc-exercises) using Visual Studio Code with the following content:

```markdown
# Bioinformatics Exercises

This repository contains my work for the Bioinformatics course, including:

- HPC cluster exercises
- Quality control
- Mapping next-generation sequences
- Variant calling
- [Additional topics to be added]

---

## 5️⃣ Prerequisites

* Completed HPC cluster exercises with local Git repository
* GitHub account (create one at github.com if needed). **The email id you use to create the account should be the same you used in git config command**

---

## 6️⃣ Exercise 4: Publishing Your Work

When you created the repository on GitHub, the default branch is main (after 2020). Look at the branch name displayed on the GitHub page. But in your local repository it is master. So lets rename the local repo branch as main before pushing the contents to github remote repository.

### Push to GitHub
```shell
# To see the default branch name
git branch
# if it is master then change it to main
git branch -M main

# push your local main branch to the remote Github repository
git push -u origin main

```

### Verify Publication
1. Visit your repository URL on GitHub
2. Verify all files are present
3. Check that .gitignore is working (no output files visible)
4. Review README formatting

---

## 7️⃣ Troubleshooting Guide

### Push Rejected
If your push is rejected due to remote changes:
```shell
git pull origin main
git push origin main
```

### Common Issues
1. **Untracked Files Appearing in Git Status**
   * Add to  .gitignore file if don't want to see that.
   * Use `git status` to verify

2. **Permission Denied**
   * Verify GitHub credentials
   * Check repository permissions

---

## 8️⃣ Overview

This goal of the exercise is to publish your course exercises to GitHub, starting with your HPC exercises and preparing for future work including mapping next-generation sequences and variant calling.

---

## 9️⃣ Environment

- All exercises are performed on the Bioinformatics server (login8.hpc.binf.unibe.ch)
- Scripts are developed and tested using SLURM job scheduler
```

### Add Documentation to Repository
```shell
git add README.md
git commit -m "Add repository documentation"
```

---

## 🔟 Exercise 5: Ongoing Workflow

### While Working on New Exercises
```shell
# Create directory for next exercise set
cd course
mkdir mapping-exercises
cd mapping-exercises

# After creating/editing files
git add .
git commit -m "Add mapping exercise solutions"
git push origin main
```

---

## 11. Best Practices

### Commit Messages
* Start with a verb (Add, Update, Fix, etc.)
* Keep first line under 50 characters
* Examples:
  * "Add variant calling scripts"
  * "Update mapping parameters for better accuracy"
  * "Fix memory allocation in SLURM script"

### Repository Organization
* Maintain separate directories for different exercise types
* Use consistent naming conventions
* Keep README.md updated
* Example structure:
  ```
  bioinformatics-exercises/
  ├── README.md
  ├── hpc-exercises/
  │   ├── scripts/
  │   └── .gitignore
  ├── mapping-exercises/
  └── variant-calling/
  ```

### Privacy and Security
* Never commit sensitive data or credentials
* Review files before committing
* Maintain an appropriate .gitignore file

---

## 12. Assessment Questions

1. What command shows configured remote repositories?
2. How do you verify .gitignore is working?
3. What is the recommended frequency for pushing changes?

---

## 13. Additional Resources

* [GitHub Documentation](https://docs.github.com)
* [Pro Git Book](https://git-scm.com/book/en/v2)
* [GitHub Guides](https://guides.github.com)
* [Ten Simple Rules for Taking Advantage of Git and GitHub](https://pmc.ncbi.nlm.nih.gov/articles/PMC4945047/#pcbi.1004947.ref005)

---

## 14. Best Practices

Before pushing:
- `git status                  # Check what's changed`
- `git pull                    # Get latest changes`
- `git add .                   # Stage changes`
- `git commit -m "message"     # Commit with clear message`
- `git push                    # Push changes`

---
